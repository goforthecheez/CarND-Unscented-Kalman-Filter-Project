#include "ukf.h"
#include "Eigen/Dense"
#include <cfloat>
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // CONSTANTS
  // px, py, v, yaw, yawd
  n_x_ = 5;
  // px, py, v, yaw, yawd, noise_a, noise_yaw
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  n_sig_ = 2 * n_aug_ + 1;
  weights_ = VectorXd::Zero(n_sig_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < n_sig_; ++i) {
    weights_(i) = 1. / (2 * (lambda_ + n_aug_));
  }

  // INSTANCE VARIABLES
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  is_initialized_ = false;
  time_us_ = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    x_ = VectorXd::Zero(n_x_);
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float r = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float rd = meas_package.raw_measurements_(2);

      x_(0) = r * cos(phi);
      x_(1) = r * sin(phi);
      float vx = rd * cos(phi);
      float vy = rd * sin(phi);
      x_(2) = sqrt(vx * vx + vy * vy);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      float px = meas_package.raw_measurements_(0);
      float py = meas_package.raw_measurements_(1);

      x_(0) = px;
      x_(1) = py;
    }

    P_ << 0.7, 0, 0, 0, 0,
          0, 0.7, 0, 0, 0,
          0, 0, 1.0, 0, 0,
          0, 0, 0, 0.3, 0,
          0, 0, 0, 0, 0.1;

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // Predict the state after motion.
  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(delta_t);

  // Update the state based on the measurement received.
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // Augment state and state covariance matrices.
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  x_aug.head(n_x_) = x_;

  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  // Generate (augmented) sigma points.
  MatrixXd A = P_aug.llt().matrixL();
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, n_sig_);
  Xsig_aug.col(0) = x_aug;
  float weight = sqrt(lambda_ + n_aug_);
  for (int i = 0; i < n_aug_; ++i) {
    Xsig_aug.col(i + 1) = x_aug + weight * A.col(i);
    Xsig_aug.col(i + n_aug_ + 1) = x_aug - weight * A.col(i);
  }

  // Pass the sigma points through the motion model, and store
  // the Xsig_pred_ values for the measurement update step.
  for (int i = 0; i < n_sig_; ++i) {
    VectorXd in = Xsig_aug.col(i);
    float p_x = in(0);
    float p_y = in(1);
    float v = in(2);
    float yaw = in(3);
    float yawd = in(4);
    float noise_a = in(5);
    float noise_yawdd = in(6);

    // Initialize out with the current values;
    VectorXd out = VectorXd(5);
    out(0) = p_x;
    out(1) = p_y;
    out(2) = v;
    out(3) = yaw;
    out(4) = yawd;

    // Set deterministic part of process model.
    if (fabs(yawd) < FLT_EPSILON) {
      out(0) += v * cos(yaw) * delta_t;
      out(1) += v * sin(yaw) * delta_t;
    } else {
      out(0) += v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      out(1) += v/yawd * (-cos(yaw + yawd * delta_t) + cos(yaw));
      out(3) += yawd * delta_t;
    }

    // Add process noise.
    float delta_t_2 = delta_t * delta_t;
    out(0) += 0.5 * delta_t_2 * cos(yaw) * noise_a;
    out(1) += 0.5 * delta_t_2 * sin(yaw) * noise_a;
    out(2) += delta_t * noise_a;
    out(3) += 0.5 * delta_t_2 * noise_yawdd;
    out(4) += delta_t * noise_yawdd;

    Xsig_pred_.col(i) = out;
  }

  // Use sigma points to estimate state mean and covariance, and use them
  // to update x_ and P_;
  x_ = VectorXd::Zero(n_x_);
  for (int i = 0; i < n_sig_; ++i) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  P_  = MatrixXd::Zero(n_x_, n_x_);
  for (int i = 0; i < n_sig_; ++i) {
    VectorXd diff = Xsig_pred_.col(i) - x_;
    // Don't forget to normalize the yaw.
    while (diff(3) > M_PI) {
      diff(3) -= 2 * M_PI;
    }
    while (diff(3) < -M_PI) {
      diff(3) += 2 * M_PI;
    }
    P_ += weights_(i) * diff * diff.transpose();
  };
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  VectorXd z = meas_package.raw_measurements_;

  // The Xsig_pred_ values have already been passed through the motion model.
  // Convert them into measurement space.
  MatrixXd Zsig = MatrixXd::Zero(2, n_sig_);
  for (int i = 0; i < n_sig_; ++i) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);

    Zsig(0, i) = px;
    Zsig(1, i) = py;
  }

  // Compute the mean predicted measurement.
  VectorXd z_pred = VectorXd::Zero(2);
  for (int i = 0; i < n_sig_; ++i) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // Compute the innovation covariance matrix S & cross-correlation
  // matrix Tc.
  MatrixXd S = MatrixXd::Zero(2, 2);
  MatrixXd Tc = MatrixXd::Zero(n_x_, 2);
  for (int i = 0; i < n_sig_; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  // Add noise to S.
  MatrixXd R = MatrixXd::Zero(2, 2);
  R(0, 0) = std_laspx_ * std_laspx_;
  R(1, 1) = std_laspy_ * std_laspy_;
  S += R;

  //Compute the Kalman gain K.
  MatrixXd K = Tc * S.inverse();

  // Update the state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  VectorXd z = meas_package.raw_measurements_;

  // The Xsig_pred_ values have already been passed through the motion model.
  // Convert them into measurement space.
  MatrixXd Zsig = MatrixXd::Zero(3, n_sig_);
  for (int i = 0; i < n_sig_; ++i) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    if (fabs(px) < FLT_EPSILON && fabs(py) < FLT_EPSILON) {
      // If px and py are both zero, phi and rho_dot will both blow up,
      // so just set this sigma point's predicted measurement to 0.
      Zsig.col(i) = VectorXd::Zero(3);
    } else {
      Zsig(0, i) = sqrt(px * px + py * py);
      Zsig(1, i) = atan2(py, px);
      Zsig(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v)/Zsig(0, i);
    }
  }

  // Compute the mean predicted measurement.
  VectorXd z_pred = VectorXd::Zero(3);
  for (int i = 0; i < n_sig_; ++i) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // Compute the innovation covariance matrix S & cross-correlation
  // matrix Tc.
  MatrixXd S = MatrixXd::Zero(3, 3);
  MatrixXd Tc = MatrixXd::Zero(n_x_, 3);
  for (int i = 0; i < n_sig_; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // Don't forget to normalize phi.
    while (z_diff(1) > M_PI) {
      z_diff(1) -= 2 * M_PI;
    }
    while (z_diff(1) < -M_PI) {
      z_diff(1) += 2 * M_PI;
    }
    S += weights_(i) * z_diff * z_diff.transpose();

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // Don't forget to normalize yaw.
    while (x_diff(3) > M_PI) {
      x_diff(3) -= 2 * M_PI;
    }
    while (x_diff(3) < -M_PI) {
      x_diff(3) += 2 * M_PI;
    }
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  // Add noise to S.
  MatrixXd R = MatrixXd::Zero(3, 3);
  R(0, 0) = std_radr_ * std_radr_;
  R(1, 1) = std_radphi_ * std_radphi_;
  R(2, 2) = std_radrd_ * std_radrd_;
  S += R;

  //Compute the Kalman gain K.
  MatrixXd K = Tc * S.inverse();

  // Update the state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  // Normalize phi.
  while (z_diff(1) > M_PI) {
    z_diff(1) -= 2 * M_PI;
  }
  while (z_diff(1) < -M_PI) {
    z_diff(1) += 2 * M_PI;
  }
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}

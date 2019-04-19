#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

// Note: This is a lightly edited version copied over from my own
// EKF project submission.
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  // Verify that a measurement was taken.
  assert(estimations.size() > 0);
  // Verify that there are as many measurements as ground truth values.
  assert(estimations.size() == ground_truth.size());

  // Compute cumulative RMSE.
  VectorXd sum = VectorXd::Zero(4);
  for (uint i = 0; i < estimations.size(); ++i) {
    VectorXd error = estimations[i] - ground_truth[i];
    VectorXd square = error.array() * error.array();
    sum += square;
  }
  VectorXd mean = sum / estimations.size();

  return sqrt(mean.array());
}

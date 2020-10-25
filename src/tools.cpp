#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  // Check the validity of inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal the ground truth vector size
  if ( (estimations.size() == 0) || (estimations.size() != ground_truth.size()) ){
    cout << "Code validity failed!" << endl;
    return rmse;
  }
  
  // Accumulate squared residuals
  for (int i = 0; i < estimations.size(); i++){
    VectorXd residual = estimations[i] - ground_truth[i];
    
    // Coefficient wise multiplication
    residual = residual.array() * residual.array();
    rmse += residual;    
  }
  
  // Calculate the mean
  rmse /= estimations.size();
  
  // Calculate the sqrt
  rmse = rmse.array().sqrt();
  
  return rmse;
  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  // check division by zero
  if ( (px == 0) && (py == 0) ){
    cout << "CalculateJacobian() - error, division by zero" << endl;
    return Hj;
  }
  
  // compute the Jacobian
  Hj(0,0) = px / sqrt( pow(px,2) + pow(py,2));
  Hj(0,1) = py / sqrt(pow(px,2) + pow(py,2));
  Hj(1,0) = -py / (pow(px,2) + pow(py,2));
  Hj(1,1) = px / (pow(px,2) + pow(py,2));
  Hj(2,0) = py * (vx*py - vy*px) / pow((pow(px,2) + pow(py,2)), 3/2);
  Hj(2,1) = px * (vy*px - vx*py) / pow((pow(px,2) + pow(py,2)), 3/2);
  Hj(2,2) = px / sqrt(pow(px,2) + pow(py,2));
  Hj(2,3) = py / sqrt(pow(px,2) + pow(py,2));
  
  return Hj;
}

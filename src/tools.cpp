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

  // Pre-compute a set of terms to avoid repeated calculation
  float c1 = px*px + py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  // Compute the Jacobian 
  Hj << (px/c2),               (py/c2),               0,     0,
       -(py/c1),               (px/c1),               0,     0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  
  return Hj;
}

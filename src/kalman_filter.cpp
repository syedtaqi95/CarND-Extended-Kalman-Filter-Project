#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft_ = F_.transpose();
  P_ = F_ * P_ * Ft_ + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd y_ = z - H_ * x_;
  MatrixXd Ht_  = H_.transpose();
  MatrixXd S_   = H_ * P_ * Ht_ + R_;
  MatrixXd Si_  = S_.inverse();
  MatrixXd K_   = P_*Ht_*Si_;

  x_ = x_ + (K_ * y_);
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K_ * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  double px = x_(0), py = x_(1), vx = x_(2), vy = x_(3);

  // Check the division by 0.
  float eps = 0.000001F;
  if (fabs(px) < eps && fabs(py) < eps) {
    px = eps;
    py = eps;
  }
  else if (fabs(px) < eps) {
    px = eps;
  }

  double rho     = sqrtf( powf(px,2) + powf(py,2) );
  double theta   = atan2(py, px);
  double rho_dot = (px * vx + py * vy) / rho;

  VectorXd z_pred(3);
  z_pred << rho, theta, rho_dot;

  VectorXd y_ = z - z_pred;

  // Keep theta within [-pi, pi]
  while(y_(1) < M_PI) y_(1) += 2*M_PI;
  while(y_(1) > M_PI) y_(1) -= 2*M_PI;

  MatrixXd Ht_  = H_.transpose();
  MatrixXd S_   = H_ * P_ * Ht_ + R_;
  MatrixXd Si_  = S_.inverse();
  MatrixXd K_   = P_ * Ht_ * Si_;

  x_ = x_ + (K_ * y_);
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K_ * H_) * P_;
}

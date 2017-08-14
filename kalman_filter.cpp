#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

	//state covariance matrix P
	P_ = MatrixXd(4, 4);
	P_ << 1, 0, 0, 0,
	      0, 1, 0, 0,
	      0, 0, 1000, 0,
	      0, 0, 0, 1000;
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;
		
  H_laser_ << 1, 0, 0, 0,
	   0, 1, 0, 0;
  F_ = MatrixXd(4, 4);
  F_ << 1, 0, 0,0,
        0, 1, 0,0,
        0, 0, 1,0,
        0,0,0,1;
}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_laser_in, MatrixXd &R_laser_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_laser_ = H_laser_in;
  R_laser_ = R_laser_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  	VectorXd z_pred = H_laser_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_laser_.transpose();
	MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_laser_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);
	double pi = 3.1415926535897;
	// Transform cartesian to polar
	float Roh_Pred=sqrt(px*px+py*py);
	
		
	float Phi_Pred=0;
        if(fabs(px)>0.0001)
	  Phi_Pred=atan2(py,px);

	float Rohdot_Pred=0;	
	if(fabs(Roh_Pred)>0.0001)			
	   Rohdot_Pred=(px*vx+py*vy)/(Roh_Pred);
	
	Hj_=tools.CalculateJacobian(x_);
	VectorXd z_pred(3);
	z_pred <<Roh_Pred,Phi_Pred,Rohdot_Pred;
	VectorXd y = z - z_pred;
	// scale phi between pi to -pi
	 
	if (y(1) > pi) 
	y(1) = fmod((y(1) - pi),(2.0*pi)) - pi;

	if (y(1) < -pi) 
	y(1) = fmod((y(1) + pi),(2.0*pi)) + pi;

	MatrixXd Ht = Hj_.transpose();
	MatrixXd S = Hj_ * P_ * Ht + R_radar_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj_) * P_;
}

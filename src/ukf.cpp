#include "ukf.h"
#include "Eigen/Dense"
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

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 30;
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
  std_yawdd_ = 1;

  
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
  
  // Set state dimension
  n_x_ = 5; 

  // Set augmented state dimension 
  n_aug_ = 7; 

  // Create vector for weights  
  weights_ = VectorXd(2*n_aug_+1); 

  // Define spreading parameters
  lambda_ = 3 - n_x_; 

  // Define sigma points matrix 
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  // matrix for sigma points in Radar measurement space
  Zsig_Radar_ = MatrixXd(3, 2 * n_aug_ + 1);


  // matrix for sigma points in Lidar measurement space
  Zsig_Lidar_ = MatrixXd(2, 2 * n_aug_ + 1);

  // Set X_sig_aug
  Xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_+1);

  // Set z_pred_ 
  z_pred_radar_ = VectorXd(3); 

  // Set S_
  S_radar_ = MatrixXd(3,3);

  // Set z_pred_ 
  z_pred_lidar_ = VectorXd(2); 

  // Set S_
  S_lidar_ = MatrixXd(2,2);

  // Set time_us_
  time_us_ = 0; 
   
}

UKF::~UKF() {}


/**
 * Process Measurements. 
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {


  if (!is_initialized_) {
  
     time_us_ = meas_package.timestamp_; 

     x_ = VectorXd(5); 
     x_ << 1, 1, 1, 1, 1; 

     P_ = MatrixXd(n_x_, n_x_); 
     P_ << 1, 0, 0, 0, 0, 
           0, 1, 0, 0, 0, 
           0, 0, 1, 0, 0, 
           0, 0, 0, 1, 0, 
           0, 0, 0, 0, 1; 
    
     if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
                
        /*
        Convert radar from spherical coordinates to cartesian coordinates and initialize state.
        I used the following web page as a guide: 
        http://tutorial.math.lamar.edu/Classes/CalcIII/SphericalCoords.aspx
        */
     
        float rho = meas_package.raw_measurements_[0]; 
        float phi = meas_package.raw_measurements_[1];     

      
        float x_cartesian = rho * cos(phi);  
        float y_cartesian = rho * sin(phi); 
      
        x_ << x_cartesian, y_cartesian, 0, 0, 0;        
      
     }


     else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        /**
        Initialize state.
        */

        float x_cartesian = meas_package.raw_measurements_[0]; 
        float y_cartesian = meas_package.raw_measurements_[1]; 
    
        x_ << x_cartesian, y_cartesian, 0, 0, 0;
     }

     
     is_initialized_ = true; 
     return;
  }

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;
  Prediction(dt);


  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
     UpdateRadar(meas_package); 
  } 

  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
     UpdateLidar(meas_package); 
  } 
    
   
}


/**
 * Generate sigma points using UKF Augmentation.
 * @param {MatrixXd*} Xsig_out in this matrix, we are going to 
 * store the augmented sigma points, we pass this matrix by reference.
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out){


   //create augmented mean vector
   VectorXd x_aug = VectorXd(n_aug_);

   //create augmented state covariance
   MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

   //create sigma point matrix
   MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

   //create augmented mean state
   x_aug.head(5) = x_;
   x_aug(5) = 0;
   x_aug(6) = 0;

   //create augmented covariance matrix
   P_aug.fill(0.0);
   P_aug.topLeftCorner(5,5) = P_;
   P_aug(5,5) = std_a_*std_a_;
   P_aug(6,6) = std_yawdd_*std_yawdd_;

   //create square root matrix
   MatrixXd L = P_aug.llt().matrixL();

   //create augmented sigma points
   Xsig_aug.col(0)  = x_aug;
   for (int i = 0; i< n_aug_; i++)
   {
      Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
   }

   *Xsig_out = Xsig_aug;

}


/**
 * Sigma Point Prediction.
 * @param {MatrixXd*} Xsig_out in this matrix, we are going to
 * store the information of the predict sigma points.
 * @param {MatrixXd} Xsig_aug the information of the augmented sigma points 
 * is contained in this matrix, we pass this matrix as a constant reference since that
 * we are not going to modify it.
 * @param {double} the delta t parameter.
 */

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, const MatrixXd &Xsig_aug, const double &delta_t){
   
   MatrixXd Xsig_pred = MatrixXd(n_x_, 2*n_aug_+1); 


   //predict sigma points
   for (int i = 0; i< 2*n_aug_+1; i++)
   {
     //extract values for better readability
     double p_x = Xsig_aug(0,i);
     double p_y = Xsig_aug(1,i);
     double v = Xsig_aug(2,i);
     double yaw = Xsig_aug(3,i);
     double yawd = Xsig_aug(4,i);
     double nu_a = Xsig_aug(5,i);
     double nu_yawdd = Xsig_aug(6,i);

     //predicted state values
     double px_p, py_p;

     //avoid division by zero
     if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
     }
     else {
        px_p = p_x + v*delta_t*cos(yaw);
       py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
   }
   
   *Xsig_out = Xsig_pred;
   
}



/**
 * Predict Mean and Covariance. 
 * @param {VectorXd*} x_out Predicted mean. 
 * @param {MatrixXd} P_out Predicted covariance.
 */

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {


  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);
  
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);



  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }

  //write result
  *x_out = x;
  *P_out = P;


}


/**
 * Predict measurement using Radar sensor.
 * @param {VectorXd*} z_out mean predicted measurement for radar. 
 * @param {MatrixXd} S_out measurement covariance matrix for radar.
 */


void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {


   //set measurement dimension, radar can measure r, phi, and r_dot
   int n_z = 3;
  
   double weight_0 = lambda_/(lambda_+n_aug_);
   weights_(0) = weight_0;
   for (int i=1; i<2*n_aug_+1; i++) {        
      double weight = 0.5/(n_aug_+lambda_);
      weights_(i) = weight;
   }
      
   //transform sigma points into measurement space
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

     // extract values for better readibility
     double p_x = Xsig_pred_(0,i);
     double p_y = Xsig_pred_(1,i);
     double v  = Xsig_pred_(2,i);
     double yaw = Xsig_pred_(3,i);

     double v1 = cos(yaw)*v;
     double v2 = sin(yaw)*v;

     // measurement model
     Zsig_Radar_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
     Zsig_Radar_(1,i) = atan2(p_y,p_x);                                 //phi
     Zsig_Radar_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
   }

   //mean predicted measurement
   VectorXd z_pred = VectorXd(n_z);
   z_pred.fill(0.0);


   for (int i=0; i < 2*n_aug_+1; i++) {
       z_pred = z_pred + weights_(i) * Zsig_Radar_.col(i);
   }

   //innovation covariance matrix S
   MatrixXd S = MatrixXd(n_z,n_z);
   S.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
     //residual
     VectorXd z_diff = Zsig_Radar_.col(i) - z_pred;

     //angle normalization
     while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
     while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

     S = S + weights_(i) * z_diff * z_diff.transpose();
   }

   //add measurement noise covariance matrix
   MatrixXd R = MatrixXd(n_z,n_z);
   R <<    std_radr_*std_radr_, 0, 0,
           0, std_radphi_*std_radphi_, 0,
           0, 0,std_radrd_*std_radrd_;

   S = S + R;

   *z_out = z_pred; 

   *S_out = S; 

}



/**
 * Predict measurement using Lidar sensor.
 * @param {VectorXd*} z_out mean predicted measurement for lidar. 
 * @param {MatrixXd} S_out measurement covariance matrix for lidar.
 */
void UKF::PredictLidarMeasurement(VectorXd* z_out, MatrixXd* S_out) {


   //set measurement dimension, lidar can measure px, py.
   int n_z = 2;
  
   double weight_0 = lambda_/(lambda_+n_aug_);
   weights_(0) = weight_0;
   for (int i=1; i<2*n_aug_+1; i++) {        
      double weight = 0.5/(n_aug_+lambda_);
      weights_(i) = weight;
   }

   //transform sigma points into measurement space
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

     // extract values for better readibility
     double p_x = Xsig_pred_(0,i);
     double p_y = Xsig_pred_(1,i);

     // measurement model
     Zsig_Lidar_(0,i) = p_x;    //p_x
     Zsig_Lidar_(1,i) = p_y;    //p_y

   }
   

   //mean predicted measurement
   VectorXd z_pred = VectorXd(n_z);
   z_pred.fill(0.0);


   for (int i=0; i < 2*n_aug_+1; i++) {
       z_pred = z_pred + weights_(i) * Zsig_Lidar_.col(i);
   }


   //innovation covariance matrix S
   MatrixXd S = MatrixXd(n_z,n_z);
   S.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
     //residual
     VectorXd z_diff = Zsig_Lidar_.col(i) - z_pred;

     //angle normalization
     while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
     while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

     S = S + weights_(i) * z_diff * z_diff.transpose();
   }


   //add measurement noise covariance matrix
   MatrixXd R = MatrixXd(n_z,n_z);
   R <<    std_laspx_*std_laspx_, 0,
           0, std_laspx_*std_laspx_;

   S = S + R;

   *z_out = z_pred; 

   *S_out = S; 
   
}


/**
 * Update state using Radar sensor.
 * @param {VectorXd} z data from Radar sensor.
 * @param {VectorXd*} x_out in this variable we are going to store the updated state mean for Radar.
 * @param {MatrixXd} P_out in this variable we are going to store the covariance matrix for Radar.
 */
void UKF::UpdateRadarState(const VectorXd &z, VectorXd* x_out, MatrixXd* P_out) {


   //set measurement dimension, radar can measure r, phi, and r_dot
   int n_z = 3;

   //set vector for weights
   double weight_0 = lambda_/(lambda_+n_aug_);
   weights_(0) = weight_0;
   for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
      double weight = 0.5/(n_aug_+lambda_);
      weights_(i) = weight;
   }
   
   //create matrix for cross correlation Tc
   MatrixXd Tc = MatrixXd(n_x_, n_z);

   //calculate cross correlation matrix
   Tc.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

     //residual
     VectorXd z_diff = Zsig_Radar_.col(i) - z_pred_radar_;
     //angle normalization
     while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
     while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

     // state difference
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     //angle normalization
     while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
     while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

     Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
   }

   //Kalman gain K;
   MatrixXd K = Tc * S_radar_.inverse();

  //residual
  VectorXd z_diff = z - z_pred_radar_;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_radar_*K.transpose();

}


/**
 * Update state using Lidar sensor.
 * @param {VectorXd} z data from Lidar sensor.
 * @param {VectorXd*} x_out in this variable we are going to store the updated state mean for Lidar.
 * @param {MatrixXd} P_out in this variable we are going to store the covariance matrix for Lidar.
 */
void UKF::UpdateLidarState(const VectorXd &z, VectorXd* x_out, MatrixXd* P_out) {


   //set measurement dimension, lidar can measure px and py
   int n_z = 2;


   //set vector for weights
   double weight_0 = lambda_/(lambda_+n_aug_);
   weights_(0) = weight_0;
   for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
      double weight = 0.5/(n_aug_+lambda_);
      weights_(i) = weight;
   }
   

   //create matrix for cross correlation Tc
   MatrixXd Tc = MatrixXd(n_x_, n_z);


   //calculate cross correlation matrix
   Tc.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

     //residual
     VectorXd z_diff = Zsig_Lidar_.col(i) - z_pred_lidar_;
     //angle normalization
     while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
     while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

     // state difference
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     //angle normalization
     while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
     while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

     Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
   }

   //Kalman gain K;
   MatrixXd K = Tc * S_lidar_.inverse();

  //residual
  VectorXd z_diff = z - z_pred_lidar_;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_lidar_*K.transpose();

}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  
  AugmentedSigmaPoints(&Xsig_aug_); 
  SigmaPointPrediction(&Xsig_pred_, Xsig_aug_, delta_t);
  PredictMeanAndCovariance(&x_, &P_); 
   
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package data from Lidar sensor.
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  PredictLidarMeasurement(&z_pred_lidar_, &S_lidar_); 
  UpdateLidarState(meas_package.raw_measurements_, &x_, &P_); 

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package data from Radar sensor.
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  PredictRadarMeasurement(&z_pred_radar_, &S_radar_); 
  UpdateRadarState(meas_package.raw_measurements_, &x_, &P_); 

}

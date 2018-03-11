#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4); 
  rmse << 0, 0, 0, 0; 
  unsigned int n = estimations.size(); 

  // We are going to check the validity of the following inputs: 
  // * the estimation vector size should not be zero
  // * the estimation vector size should equal ground truth size 

  if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
     
     cout << "Invalid estimation or ground_thruth data" << endl; 
     return rmse; 

  }

  // accumulate square residuals 
  for (unsigned int i=0; i < estimations.size(); ++i)
  {
     VectorXd residual = estimations[i] - ground_truth[i]; 
     
     residual = residual.array() * residual.array(); 
     rmse += residual; 

  }


  // calculate the mean
  rmse = rmse / n; 

  // calculate the squared root 
  rmse = rmse.array().sqrt(); 

  // return the result 
  return rmse; 


}

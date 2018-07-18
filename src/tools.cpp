#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd residuals(4);
	residuals << 0,0,0,0;


	if (estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return residuals;
	}
    
	for(int i=0; i < estimations.size(); ++i){
        VectorXd diff = estimations[i] - ground_truth[i];
		residuals += VectorXd(diff.array()*diff.array());
	}

	residuals /= estimations.size();
	residuals = residuals.array().sqrt();

	return residuals;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//check division by zero
	if (px + py == 0) 
	    cout << "Calculate Jacobion - Error - Division by Zero";
	
	//compute the Jacobian matrix
	float px2 = px*px;
	float py2 = py*py;
	Hj << px/sqrt(px2+py2), py/sqrt(px2+py2), 0, 0,
	      -py/(px2+py2), px/(px2+py2), 0, 0,
	      py*(vx*py-vy*px)/pow((px2+py2),3/2), px*(vy*px-vx*py)/pow((px2+py2),3/2), px/sqrt(px2+py2), py/sqrt(px2+py2);

	return Hj;
}

/*
Filename: solve.cpp
Description: Solve the deterministic or stochastic ODE 
			 for the gene regulatory network.
Author: Lei Zhao	Version: 2.0	Date: 2016-9-20
*/

#include "../include/solve.h"


/*
Euler method to solve the deterministic differential equations.
:param dx: the differential value of x with the force function.
:param x: a set of values corresponding to the expression of every gene in model.
:param force: the instance of class Force, include the instance of class Net.
:return: false if the result x value out of the boundary.
*/
bool euler_det(Vec_DP &dx, Vec_DP &x, Force &force) {
    int i;
    force.cal_dx(dx, x);
	// update the x value
    for (i=0; i<x.size(); i++) {
    	x[i] += EULER_STEP * dx[i];
    }
    if (is_out_bound(x, force.x_min, force.x_max)) {
    	//cout << "The x value is out of boundary!!!" << endl;
    	return false;
    }
    return true;
}


/*
Euler method with noise term to solve the stochastic differential equations.
:param dx: the differential value of x with the force function.
:param x: a set of values corresponding to the expression of every gene in model.
:param force: the instance of class Force, include the instance of class Net.
:return: false if the result x value out of the boundary, the x do not change if return false.
*/
bool euler_sto(Vec_DP &dx, Vec_DP &x, Force &force) {
	int i;
	extern int idum;
	Vec_DP x_tmp = x;
	force.cal_dx(dx, x);
	for (i=0; i<x.size(); i++) {
    	x_tmp[i] += EULER_STEP * dx[i] + sqrt(2*DIFF*EULER_STEP)*gasdev(idum);
    	// check the boundry conditions
    	if (is_out_bound(x_tmp, force.x_min, force.x_max)) {
    		//cout << "The x value is out of boundary!!!" << endl;
    		return false;
    	}
	}
	x = x_tmp;		
	return true;
}



/*
Filename: force.h
Description: Define the force field that used to model the network dynamics.
Author: Lei Zhao	Version: 2.0	Date: 2016-9-20
*/


#ifndef _FORCE_H_
#define _FORCE_H_

#include <cmath>
#include "../include/net.h"

typedef Vec<double> Vec_DP;
typedef Mat<double> Mat_DP;
typedef Vec<int> Vec_INT;
typedef Mat<int> Mat_INT;


class Force {
public:
	// default constructor
	Force(const Net &net);
	// init the parameters used in Hill equation
	void init_hill_para();
	// init the parameters used in ode1_force
	void init_ode1_para();
	// destruct function
	~Force();
	
	
	/*
	Wrapper of force function, which provide a unified calling interface
	:para dx: the differential value of x with the force function.
	:param x: a set of values corresponding to the expression of every gene in model.
	*/
	void cal_dx(Vec_DP &dx, const Vec_DP &x);
	
	/*
	Force function ODE to calculate the dynamics of the model, every regulation is
	Hill method formatted. Use the additive rule mentioned in (Chao Tang, Quantitative Biology, 2013)
	:para dx: the differential value of x with the force function.
	:param x: a set of values corresponding to the expression of every gene in model.
	*/
	void cal_dx_hill(Vec_DP &dx, const Vec_DP &x);
	
	/*
	Force function ode1_force is mentioned in paper(Tyson and Novak, 2010)
	:para dx: the differential value of x with the force function.
	:param x: a set of values corresponding to the expression of every gene in model.
	*/
	void cal_dx_ode1(Vec_DP &dx, const Vec_DP &x);
	
	
	// instance of Net class
	Net net;
	
	// the boundary of x value
	double x_min;
	double x_max;
	
	// parameters for ode1 force
	// steepness of the sigmoidal function, 1<=sigma<=20
	Vec_DP sigma_vec;
	// timescale, 0.1<=gamma<=10, default=1
	Vec_DP gamma_vec;
	// leak weight, default=0
	Vec_DP w_leak_vec;
	
	// parameters for Hill force 
	// Hill coefficient constant
	double n;
	// degradation rate
	double k;
	// regulation strength
	Mat_DP s_mat;
	// the x power n Hill coefficient
	Vec_DP xn;
	// the s power n Hill coefficient
	Mat_DP sn_mat;

	
};




#endif

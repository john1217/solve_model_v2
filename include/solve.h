/*
Filename: solve.h
Description: Solve the deterministic or stochastic ODE 
			 for the gene regulatory network.
Author: Lei Zhao	Version: 2.0	Date: 2016-9-20
*/

#ifndef _SOLVE_H_
#define _SOLVE_H_

#include <iostream>
#include "../include/basic.h"
#include "../include/force.h"
#include "../include/rand_num.h"

using namespace std;

typedef Vec<double> Vec_DP;
typedef Mat<double> Mat_DP;
typedef Vec<int> Vec_INT;
typedef Mat<int> Mat_INT;

// the step length for euler method
const double EULER_STEP = 0.01;
//diffusion coefficient
const double DIFF = 0.01;


/*
Euler method to solve the deterministic differential equations.
:param dx: the differential value of x with the force function.
:param x: a set of values corresponding to the expression of every gene in model.
:param force: the instance of class Force, include the instance of class Net.
:return: false if the result x value out of the boundary.
*/
bool euler_det(Vec_DP &dx, Vec_DP &x, Force &force);


/*
Euler method with noise term to solve the stochastic differential equations.
:param dx: the differential value of x with the force function.
:param x: a set of values corresponding to the expression of every gene in model.
:param force: the instance of class Force, include the instance of class Net.
:return: false if the result x value out of the boundary, the x do not change if return false.
*/
bool euler_sto(Vec_DP &dx, Vec_DP &x, Force &force);




/*
	// parameters for calculate fixed points
	// number of initial value
	int ninit;
	// max number of fixed point 
	int nfp_max;
	// parameters for calculate landscape
	// the number of grid 
	int ngrid;
	// the max number of step to calcualte
	int nstep_max;
	// map landscape into two dimensions
	int dim_x;
	int dim_y;
	
	
	// map landscape into two dimensions
	dim_x = 0;	//DAF-16
	dim_y = 1;	//TORC1
	// parameters for calculate fixed points
	ninit = (int)1e3;
	// largest number of fixed points
	nfp_max = 5;
	// the number of grid 
	ngrid = 200;
	// the max number of step to calcualte
	nstep_max = (int)1e8;
*/

#endif



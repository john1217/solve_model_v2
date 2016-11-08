/*
Filename: force.cpp
Description: Define the force field that used to model the network dynamics.
Author: Lei Zhao	Version: 2.0	Date: 2016-9-20
*/

#include "../include/force.h"


// default constructor
Force::Force(const Net &net) {
	this->net = net;
	init_hill_para();
	init_ode1_para();
}

// init the parameters used in ode1_force
void Force::init_ode1_para() {
	sigma_vec = Vec_DP(1.0, net.n_node);
	gamma_vec = Vec_DP(1.0, net.n_node);
	w_leak_vec = Vec_DP(0.0, net.n_node);
}

// init the parameters used in Hill equation
void Force::init_hill_para() {
	// Hill coefficient
	n = 3.0;
	// degradation rate
	k = 1.0;
	// the boundary of x value
	x_min = 0.0;
	x_max = 1.0;
	// the x power n Hill coefficient
	xn = Vec_DP(0.0, net.n_node);
	// regulation strength
	s_mat = Mat_DP(0.0, net.n_node, net.n_node);
	// the s power n Hill coefficient
	sn_mat = Mat_DP(0.0, net.n_node, net.n_node);
	for (int i=0; i<net.n_node; i++) {
		for (int j=0; j<net.n_node; j++) {
			if (net.t_mat[i][j] != 0) {
				s_mat[i][j] = 0.5;
				sn_mat[i][j] = pow(s_mat[i][j],n);
			}
		}
	}
}


/*
Wrapper of calculate dx function, which provide a unified calling interface
:para dx: the differential value of x with the force function.
:param x: a set of values corresponding to the expression of every gene in model.
*/
void Force::cal_dx(Vec_DP &dx, const Vec_DP &x) {
	cal_dx_ode1(dx, x);
	//cal_dx_hill(dx, x);
}


/*
Force function cal_dx_ode1 is mentioned in paper(Tyson and Novak, 2010)
:para dx: the differential value of x with the force function.
:param x: a set of values corresponding to the expression of every gene in model.
*/
void Force::cal_dx_ode1(Vec_DP &dx, const Vec_DP &x) {
	int i, j;
	double w;
	dx = 0.0;
	for (j=0; j<net.n_node; j++) {
		w = 0.0;
		for (i=0; i<net.n_node; i++) {
			w += net.t_mat[i][j] * net.w_mat[i][j] * x[i];
		}
		w += w_leak_vec[j];
		dx[j] = gamma_vec[j] * (1.0/(1.0+exp(-sigma_vec[j]*w))-x[j]);
	}
}


/*
Force function ODE to calculate the dynamics of the model, every regulation is
Hill method formatted. Use the additive rule mentioned in (Chao Tang, Quantitative Biology, 2013)
:para dx: the differential value of x with the force function.
:param x: a set of values corresponding to the expression of every gene in model.
*/
void Force::cal_dx_hill(Vec_DP &dx, const Vec_DP &x) {	
	int i, j;
	dx = 0.0;
	for (i=0; i<net.n_node; i++) {
		xn[i] = pow(x[i], n);
	}
	for (j=0; j<net.n_node; j++) {
		for (i=0; i<net.n_node; i++) {
			if (net.t_mat[i][j] > 0) {
				dx[j] += net.w_mat[i][j] * xn[i]/(sn_mat[i][j]+xn[i]);
			} else if (net.t_mat[i][j] < 0) {
				dx[j] += net.w_mat[i][j] * sn_mat[i][j]/(sn_mat[i][j]+xn[i]);
			}
		}
		dx[j] = dx[j] - k * x[j];
	}
}

//destruct function
Force::~Force() {
}

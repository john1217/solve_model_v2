/*
Filename: analysis.cpp
Description: Analysis the dynamics of the network model.
Author: Lei Zhao	Version: 2.0	Date: 2016-9-21
*/

#include "../include/analysis.h"

/*
Find a solution of given x for the deterministic ODE.
:param x: the initial values corresponding to every gene in model.
:param force: the instance of class Force.
*/
void find_solution(Vec_DP &x, Force &force) {
	Vec_DP dx = Vec_DP(0.0, x.size());
	for (int i=0; i<MAX_STEP; i++) {
		euler_det(dx, x, force);
		if (is_zero(dx, 1e-6)) {
			return;
		}
	}
	cout << "Cannot reach the fixed point in MAX_STEP!!!" << endl;
}
	
/*
Find all solutions of the deterministic ODE.
:param solutions: get all the solutions and saved in this matrix.
:param force: the instance of class Force include the system ODE.
:param n_init: the number of initial values to find the solutions
*/
void find_all_solutions(Mat_DP &solutions, Force &force, const int n_init) {
	Vec_DP x = Vec_DP(0.0, force.net.n_node);
	Vec_DP dx = x;
	Mat_DP all_x = Mat_DP(0.0, n_init, x.size());
	for (int i=0; i<n_init; i++) {
		rand_vec(x, force.x_min, force.x_max);
		find_solution(x, force);
		for (int j=0; j<x.size(); j++) {
			all_x[i][j] = x[j];
		}
	}
	unique_mat_row(solutions, all_x);
}

/*
Create the frequency landscape by putting N trajectory component into land.
:param land: the frequency landscape, the bin edges along the x and y dimension
:param fx, fy: the integral of forces in the 2D landscape for column xid, yid 
:param xid, yid: the column ids used to generate the landscape
:param n_step: the number of steps we want to map to the 2D landscape,
:param force: the instance of class Force include the system ODE.
*/
void freq_land(Mat_INT &land, Mat_DP &fx, Mat_DP &fy, int xid, int yid, int n_step, Force &force) {
	int i, xx, yy;
	int n_node = force.net.n_node;
	int n_grid = land.nrows();
	double min = force.x_min;
	double max = force.x_max;
	int elaspe = 1000;
	// the interval between each grid
	double delta = (max-min)/(n_grid-1);
	Vec_DP dx(0.0,n_node);
	// initial x value
	Vec_DP x_value(0.0,n_node);
	rand_vec(x_value, min, max);
	// make the x_value stable
	for (i=0; i<elaspe; i++) {
		euler_det(dx, x_value, force);
	}
	land = 0;
	fx = 0;
	fy = 0;
	for (i=0; i<n_step; i++) {
		if (euler_sto(dx, x_value, force)) {
			xx = (int)((x_value[xid]-min)/delta);
			yy = (int)((x_value[yid]-min)/delta);
			land[yy][xx]++;
			fx[yy][xx] += dx[xid];
			fy[yy][xx] += dx[yid];
		}
	}
}

















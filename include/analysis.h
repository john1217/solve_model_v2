/*
Filename: analysis.h
Description: Analysis the dynamics of the network model.
Author: Lei Zhao	Version: 2.0	Date: 2016-9-21
*/

#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

#include "../include/force.h"
#include "../include/basic.h"
#include "../include/solve.h"

const int MAX_STEP = (int)1e8;

/*
Find a solution of given x for the deterministic ODE.
:param x: the initial values corresponding to every gene in model.
:param force: the instance of class Force.
*/
void find_solution(Vec_DP &x, Force &force);

/*
Find all solutions of the deterministic ODE.
:param solutions: get all the solutions and saved in this matrix.
:param force: the instance of class Force include the system ODE.
:param n_init: the number of initial values to find the solutions
*/
void find_all_solutions(Mat_DP &solutions, Force &force, const int n_init);

/*
Create the frequency landscape by putting N trajectory component into land.
:param land: the frequency landscape, the bin edges along the x and y dimension
:param fx, fy: the integral of forces in the 2D landscape for column xid, yid 
:param xid, yid: the column ids used to generate the landscape
:param n_step: the number of steps we want to map to the 2D landscape,
:param force: the instance of class Force include the system ODE.
*/
void freq_land(Mat_INT &land, Mat_DP &fx, Mat_DP &fy, int xid, int yid, int n_step, Force &force);

#endif

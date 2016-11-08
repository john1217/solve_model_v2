/*
Filename: traj.cpp
Description: Generate and manipulate the trajectory
Author: Lei Zhao	Version: 2.0	Date: 2016-9-20
*/

#include "../include/trajectory.h"


/*
Generate a deterministic trajectory for the given size and start point
:param traj: the matrix with a given size to store the trajectory
:param x_start: the start point of the trajectory
:param column_ids: the variable ids to store in traj, the length of column_ids should be equal to the ncol of traj
:param force: the instance of class Force used to solve the equations
:return: if return smaller than ncol of matrix, the return means the steps used to reach the fixed point
*/
int generate_traj_det(Mat_DP &traj, const Vec_DP x_start, const Vec_INT &column_ids, Force &force) {
	int i, j;
	double tmp;
	int len = traj.nrows();
	Vec_DP x = x_start;
	Vec_DP dx = Vec_DP(0.0, x.size());
	// give the first row of trajectory as the x_start value
	for (j=0; j<column_ids.size(); j++) {
		traj[0][j] = x_start[column_ids[j]];
	}
	for (i=1; i<len; i++) {
		if (euler_det(dx, x, force)) {
			for (j=0; j<column_ids.size(); j++) {
				traj[i][j] = x[column_ids[j]];
			}
			// if dx=0, means stable, return the step id when the system reach the fixed point
			if (is_zero(dx, 1e-6)) {
				return i;
			}
		}
	}
	return i;
}


/*
Generate a stochastic trajectory for the given size and start point
:param traj: the matrix with a given size to store the trajectory
:param x_start: the start point of the trajectory
:param column_ids: the variable ids to store in traj, the length of column_ids should be equal to the ncol of traj
:param force: the instance of class Force used to solve the equations
:return: the length of the trajectory matrix
*/
int generate_traj_sto(Mat_DP &traj, const Vec_DP x_start, const Vec_INT &column_ids, Force &force) {
	int i, j;
	int len = traj.nrows();
	Vec_DP x = x_start;
	Vec_DP dx = Vec_DP(0.0, x.size());
	// give the first row of trajectory as the x_start value
	for (j=0; j<column_ids.size(); j++) {
		traj[0][j] = x_start[column_ids[j]];
	}
	for (i=1; i<len; i++) {
		while (!euler_sto(dx, x, force));
		for (j=0; j<column_ids.size(); j++) {
			traj[i][j] = x[column_ids[j]];
		}
	}
	return i;
}





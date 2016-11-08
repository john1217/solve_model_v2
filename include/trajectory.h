/*
Filename: traj.h
Description: Generate and manipulate the trajectory
Author: Lei Zhao	Version: 2.0	Date: 2016-9-20
*/


#ifndef _TRAJECTORY_H_
#define _TRAJECTORY_H_

#include "../include/solve.h"
#include "../include/basic.h"
#include "../include/force.h"

using namespace std;

typedef Vec<double> Vec_DP;
typedef Mat<double> Mat_DP;
typedef Vec<int> Vec_INT;
typedef Mat<int> Mat_INT;



/*
Generate a deterministic trajectory for the given size and start point
:param traj: the matrix with a given size to store the trajectory
:param x_start: the start point of the trajectory
:param column_ids: the variable ids to store in traj, the length of column_ids should be equal to the ncol of traj
:param force: the instance of class Force used to solve the equations
:return: if return smaller than ncol of matrix, the return means the steps used to reach the fixed point
*/
int generate_traj_det(Mat_DP &traj, const Vec_DP x_start, const Vec_INT &column_ids, Force &force);


/*
Generate a stochastic trajectory for the given size and start point
:param traj: the matrix with a given size to store the trajectory
:param x_start: the start point of the trajectory
:param column_ids: the variable ids to store in traj, the length of column_ids should be equal to the ncol of traj
:param force: the instance of class Force used to solve the equations
:return: the length of the trajectory matrix
*/
int generate_traj_sto(Mat_DP &traj, const Vec_DP x_start, const Vec_INT &column_ids, Force &force);



#endif

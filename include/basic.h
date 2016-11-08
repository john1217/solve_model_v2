/*
Filename: basic.h
Description: Basic tools used in the calculation
Author: Lei Zhao	Version: 2.0	Date: 2016-6-23
*/

#ifndef _BASIC_H_
#define _BASIC_H_

#include <ctime>
#include <cstring>
//#include <io.h>
//#include <direct.h>
#include <unistd.h>
#include <fcntl.h>

#include "../include/data_struct.h"
#include "../include/net.h"
#include "../include/rand_num.h"

typedef Vec<double> Vec_DP;
typedef Mat<double> Mat_DP;
typedef Vec<int> Vec_INT;
typedef Mat<int> Mat_INT;


/*
Delete the redundant row of the given matrix
:param unique_mat: initialize a new matrix which contains the unique row values
:param mat: the input matrix to delete redundant row.
*/
void unique_mat_row(Mat_DP &unique_mat, const Mat_DP &mat);

// random the x value in the boundary condition
void rand_vec(Vec_DP &x, const double x_min, const double x_max);

// Return true if any element in x cross the boundary
bool is_out_bound(Vec_DP &x, const double x_min, const double x_max);

// Return true if all element in x is zero
bool is_zero(Vec_DP &x, const double tolerance);

// Return true if all element in x1 is equal to x2
bool is_equal(Vec_DP &x1, Vec_DP &x2, const double tolerance);

// Return true if the id1 row in mat1 is equal to the id2 row in mat2
bool is_equal_mat_row(Mat_DP &mat1, const int id1, Mat_DP &mat2, const int id2, const double tolerance);

// transform the vector to matrix by the topology of the network
void vec2mat(Mat_DP &m, Vec_DP &v, const Net &net);

// transform the matrix to vector by the topology of the network	
void mat2vec(Vec_DP &v, Mat_DP &m, const Net &net);

// print the current time
void printTime();

// print the current day
string printDay();

// make dir if the dir is not exist
void makeDir(string fname);

#endif

/*
Filename: basic.cpp
Description: Basic tools used in the calculation
Author: Lei Zhao	Version: 2.0	Date: 2016-6-23
*/

#include "../include/basic.h"

/*
Delete the redundant row of the given matrix
:param unique_mat: initialize a new matrix which contains the unique row values
:param mat: the input matrix to delete redundant row.
*/
void unique_mat_row(Mat_DP &unique_mat, const Mat_DP &mat) {
	int i, j, k;
	// tolerance to compare the value
	double const tolerance = 0.01;
	Mat_DP mat_tmp = mat;
	int n_unique = 1;
	int n_col = mat.ncols();
	for (i=1; i<mat.nrows(); i++) {
		for (j=0; j<mat.nrows(); j++) {
			if (!is_equal_mat_row(mat_tmp, j, mat_tmp, i, tolerance)) {
				if (j < n_unique-1) {
					continue;
				} else {
					for (k=0; k<n_col; k++) {
						mat_tmp[n_unique][k] = mat_tmp[i][k];
					}
					n_unique++;
					break;
				}
			} else {
				break;
			}
		}
	}
	unique_mat = Mat_DP(0.0, n_unique, n_col);
	for (i=0; i<n_unique; i++) {
		for (j=0; j<n_col; j++) {
			unique_mat[i][j] = mat_tmp[i][j];
		}
	}
}

// random the x value in the boundary condition
void rand_vec(Vec_DP &x, const double x_min, const double x_max ) {
	extern int idum;
	for (int i=0; i<x.size(); i++) {
		x[i] = ran1(idum)*(x_max-x_min) + x_min;
	}
}

// Return true if any element in x cross the boundary
bool is_out_bound(Vec_DP &x, const double x_min, const double x_max) {
	// check every element of x
	for (int i=0; i<x.size(); i++) {
		if (x[i]<x_min || x[i]>x_max) {
			return true;
		}
	}
	return false;
}

// Return true if all element in x is zero
bool is_zero(Vec_DP &x, const double tolerance) {
	// the tolerance to test if x is zero
	for (int i=0; i<x.size(); i++) {
		if (x[i] < -tolerance || x[i] > tolerance) {
		return false;
		}
	}
	return true;
}

// Return true if all element in x1 is equal to x2
bool is_equal(Vec_DP &x1, Vec_DP &x2, const double tolerance) {
	double dx;
	for (int i=0; i<x1.size(); i++) {
		dx = x1[i] - x2[i]; 
		if (dx < -tolerance || dx > tolerance) {
			return false;
		}
	}
	return true;
}

// Return true if the id1 row in mat1 is equal to the id2 row in mat2
bool is_equal_mat_row(Mat_DP &mat1, const int id1, Mat_DP &mat2, const int id2, const double tolerance) {
	double dx;
	for (int i=0; i<mat1.ncols(); i++) {
		dx = mat1[id1][i] - mat2[id2][i];
		if (dx < -tolerance || dx > tolerance) {
			return false;
		}
	}
	return true;
}

// transform the vector to matrix by the topology of the network	
void vec2mat(Mat_DP &m, Vec_DP &v, const Net &net){
	int dim = m.nrows();
	int k = 0;
	for (int i=0; i<dim; i++) {
		for (int j=0; j<dim; j++) {
			if (net.t_mat[i][j] != 0) {
				m[i][j] = v[k];
				k++;
			}
		}
	}
}

// transform the matrix to vector by the topology of the network	
void mat2vec(Vec_DP &v, Mat_DP &m, const Net &net){
	int dim = m.nrows();
	int k = 0;
	for (int i=0; i<dim; i++) {
		for (int j=0; j<dim; j++) {
			if (net.t_mat[i][j] != 0) {
				v[k] = m[i][j]; 
				k++;
			}
		}
	}
}

// print the current time
void printTime() {
	time_t now = time(NULL);
	tm *ptm = localtime(&now);
	char ctm[26];
	strftime(ctm,26,"%Y-%m-%d %H:%M:%S",ptm);
	cout << ctm << endl;
}

// print the current day
string printDay() {
	time_t now = time(NULL);
	tm *ptm = localtime(&now);
	char ctm[26];
	strftime(ctm,26,"%Y-%m-%d",ptm);
	return string(ctm);
}

// make dir if the dir is not exist
void makeDir(string fname) {
	char filename[1000];
	strcpy(filename,fname.c_str());
	char *tag;
	for (tag=filename; *tag; tag++) {
		if (*tag == '/') {
			char buf[1000];
			char path[1000];
			strcpy(buf,filename);
			buf[strlen(filename)-strlen(tag)+1] = (char)NULL;
			strcpy(path,buf);
			if (access(path,6) == -1) {
				mkdir(path);
			}
		}
	}
}





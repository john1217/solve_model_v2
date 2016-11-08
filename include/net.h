/*
Filename: net.h
Description: Define the topology and related parameters like the weight of 
			 each regulations for the gene regulatory network.
Author: Lei Zhao	Version: 2.0	Date: 2016-9-20
*/


#ifndef _NET_H_
#define _NET_H_


#include "../include/data_struct.h"
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>

using namespace std;

typedef Vec<double> Vec_DP;
typedef Mat<double> Mat_DP;
typedef Vec<int> Vec_INT;
typedef Mat<int> Mat_INT;

class Net {
public:
	// default constructor
	Net();
	// construction from file matrix
	Net(string f_t_mat, string f_gene_symbols, int n_node);
	// copy constructor
	Net(const Net &net);
	// scale the w_mat, for every gene, the sum of regulation weight equal to 1
	void scale_weight();
	// init the run parameters
	void init_para();
	// destruct function
	~Net();
	
	// name of every factor 
	string *gene_symbols;
	// model topology matrix
	Mat_INT t_mat;
	// model weight matrix, need to be scaled
	Mat_DP w_mat;
	// network dim, number of nodes
	int n_node;
	// number of edges in network
	//int n_edge;
	// in degree for each nodes
	//Vec_INT indegree;
};

#endif

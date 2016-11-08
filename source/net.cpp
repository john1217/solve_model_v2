/*
Filename: net.cpp
Description: Define the topology and related parameters like the weight of 
			 each regulations for the gene regulatory network.
Author: Lei Zhao	Version: 2.0	Date: 2016-9-20
*/


#include "../include/net.h"


// default constructor
Net::Net() {
}

// construction from file matrix
Net::Net(string f_t_mat, string f_gene_symbols, int n_node) {
	this->n_node = n_node;
	ifstream infile_mat(f_t_mat.c_str());
	// divide the matrix into adjacency and strength matrix
	t_mat = Mat_INT(infile_mat,n_node,n_node);
	infile_mat.close();
	w_mat = Mat_DP(0.0,n_node,n_node);
	for (int i=0; i<n_node; i++) {
		for (int j=0; j<n_node; j++) {
			if (t_mat[i][j] != 0) {
				//s_mat[i][j] = 0.5;
				w_mat[i][j] = 1.0;
			}
		}
	}
	//scale the weight matrix
	//scale_weight();
	gene_symbols = new string[n_node];
	ifstream infile_node(f_gene_symbols.c_str());
	for (int i=0; i<n_node; i++) {
		infile_node >> gene_symbols[i];
	}
	infile_node.close();
	init_para();
	//cout << "Net init OK!" << endl;
}

// copy constructor
Net::Net(const Net &net) {
	n_node = net.n_node;
	//n_edge = net.n_edge;
	gene_symbols = new string[n_node];
	for (int i=0; i<n_node; i++) {
		gene_symbols[i] = net.gene_symbols[i];
	}
	t_mat = net.t_mat;
	w_mat = net.w_mat;
	init_para();
}

// scale the w_mat, for every gene, the sum of regulation weight equal to 1
void Net::scale_weight() {
	Vec_DP scale(0.0, n_node);
	for (int i=0; i<n_node; i++) {
		for (int j=0; j<n_node; j++) {
			// regulation from i to j
			if (abs(w_mat[i][j]) > 1e-6) {
				scale[j] += abs(w_mat[i][j]);
			}
		}
	}
	for (int i=0; i<n_node; i++) {
		for (int j=0; j<n_node; j++) {
			if (abs(scale[j]) > 1e-6) {
				w_mat[i][j] = w_mat[i][j]/scale[j];
			} else {
				w_mat[i][j] = 0;
			}
		}
	}
} 

// init the run parameters
void Net::init_para() {
	/*
	// number of edges
	n_edge = 0;
	//in degree for each node
	indegree = Vec_INT(0, n_node);
	for (int i=0; i<n_node; i++) {
		for (int j=0; j<n_node; j++) {
			if (abs(t_mat[i][j]) > 1e-6) {
				n_edge++;
				indegree[j]++;
			}
		}
	}
	*/
}

// destruct function
Net::~Net() {
	//delete[] gene_symbols;
}



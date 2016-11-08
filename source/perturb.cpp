/*
Filename: perturb.cpp
Description: Perturb the parameters of the network or force 
Author: Lei Zhao	Version: 2.0	Date: 2016-10-1
*/

#include "../include/perturb.h"


/*
Remove the selected edge of the network.
:param net: the instance of class Net.
:param src_id: the source node id of the regulation need to remove.
:param tag_id: the target node id of the regulation need to remove.
*/
void remove_edge(Net &net, const int src_id, const int tag_id) {
	net.t_mat[src_id][tag_id] = 0;
	net.w_mat[src_id][tag_id] = 0.0;
	//scale the weight matrix
	//net.scale_weight();
}

/*
Random given the value of the weight matrix
:param net: the instance of class Net.
:param w_min: the min value of random weight.
:param w_max: the max value of random weight.
*/
void rand_w_mat(Net &net, const double w_min, const double w_max) {
	extern int idum;
	for (int i=0; i<net.n_node; i++) {
		for (int j=0; j<net.n_node; j++) {
			if (net.w_mat[i][j] != 0.0) {
				net.w_mat[i][j] = ran1(idum)*(w_max-w_min) + w_min;
			}
		}
	}
}

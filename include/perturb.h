/*
Filename: perturb.h
Description: Perturb the parameters of the network or force 
Author: Lei Zhao	Version: 2.0	Date: 2016-10-1
*/

#ifndef _PERTURB_H_
#define _PERTURB_H_


#include "../include/net.h"
#include "../include/rand_num.h"

/*
Remove the selected edge of the network.
:param src_id: the source node id of the regulation need to remove.
:param tag_id: the target node id of the regulation need to remove.
:param net: the instance of class Net.
*/
void remove_edge(Net &net, const int src_id, const int tag_id);


/*
Random given the value of the weight matrix
:param net: the instance of class Net.
:param w_min: the min value of random weight.
:param w_max: the max value of random weight.
*/
void rand_w_mat(Net &net, const double w_min, const double w_max);


#endif

/*
Filename: rand_num.h
Description: random number generation 
Author: Lei Zhao	Version: 0.1	Date: 2014-04-07	
*/
#ifndef _RAND_NUM_H_
#define	_RAND_NUM_H_

#include <ctime>
#include <cstdlib>
#include <cmath>


using namespace std;

// random number in (0,1) generation from "Numberical Recipes in C++" P208 
double ran1(int &idum);

// gauss random number, mean=0, sd=1, generation "Numberical Recipes in C++" P216
double gasdev(int &idum);

#endif

/*
Filename: rand_num.cpp
Description: random number generation 
Author: Lei Zhao	Version: 0.1	Date: 2014-04-07	
*/
#include "../include/rand_num.h"

// random number in (0,1) get from "Numberical Recipes in C++" P208 
double ran1(int &idum) {
	const int IA=16807, IM=2147483647, IQ=127773, IR=2836, NTAB=32;
	const int NDIV = (1+(IM-1)/NTAB);
	const double EPS=3.0e-16, AM=1.0/IM, RNMX=(1.0-EPS);
	static int iy=0;
	//static Vec_INT iv(NTAB);
	static int iv[NTAB];
	int j, k;
	double temp;
	if (idum <= 0 || !iy) {
		if (-idum < 1) idum=1;
		else idum = -idum;
		for (j=NTAB+7;j>=0;j--) {
			k=idum/IQ;
			idum=IA*(idum-k*IQ)-IR*k;
			if (idum < 0) idum += IM;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}


// gauss random number mean=0, sd=1, get from "Numberical Recipes in C++" P216
double gasdev(int &idum) {
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;
	if (idum < 0) iset=0;
	if (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

// randomly generate 0 or 1, from "Numberical Recipes in C++" P221
int irbit1(unsigned long &iseed) {
	unsigned long newbit;
	newbit = ((iseed>>17)&1)^((iseed>>4)&1)^((iseed>>1)&1)^(iseed&1);
	iseed = (iseed << 1) | newbit;
	return int(newbit);
}

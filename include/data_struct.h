/*
Filename: data_struct.h
Description: data structure and related functions about vector and matrix derived from "Numerical Recipes in C++"
Author: Lei Zhao	Version: 0.2	Date: 2014-04-10	
*/

#ifndef _DATA_STRUCT_H_
#define _DATA_STRUCT_H_

#include <string>
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
using namespace std;

template <class T>
class Mat;
template <class T>
class Vec;

template <class T>
class Vec {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	Vec();
	explicit Vec(int n);		// Zero-based array
	Vec(const T &a, int n);	//initialize to constant value
	Vec(const T *a, int n);	// initialize to array
	Vec(ifstream &f, int n);	//initialize to file
	Vec(const Vec &rhs);	// Copy constructor
	Vec & operator=(const Vec &rhs);	//assignment
	Vec & operator=(const T &a);	//assign a to every element
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	T min() const;
	T max() const;
	T dot_product(const Vec &rhs) const;	//dot product to rhs
	double dist(const Vec &rhs) const;	// distance between the two vectors
	void print() const;
	void write(ofstream &f) const;
	void write(string fname) const;
	// initialize from the id th row(dim:1) or col(dim:2) of a matrix 
	void from_mat(const Mat<T> &mat, int dim, int id); 
	double mean() const;
	double var() const;		//variance
	~Vec();
};
template <class T>
Vec<T>::Vec() : nn(0), v(0) {}

template <class T>
Vec<T>::Vec(int n) : nn(n), v(new T[n]) {}

template <class T>
Vec<T>::Vec(const T& a, int n) : nn(n), v(new T[n])
{
	for(int i=0; i<n; i++)
		v[i] = a;
}

template <class T>
Vec<T>::Vec(const T *a, int n) : nn(n), v(new T[n])
{
	for(int i=0; i<n; i++)
		v[i] = *a++;
}

template <class T>
Vec<T>::Vec(ifstream &f, int n) : nn(n), v(new T[n]) 
{
	for (int i=0; i<n; i++)
		f >> v[i];
}


template <class T>
Vec<T>::Vec(const Vec<T> &rhs) : nn(rhs.nn), v(new T[nn])
{
	for(int i=0; i<nn; i++)
		v[i] = rhs[i];
}

template <class T>
Vec<T> & Vec<T>::operator=(const Vec<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (nn != rhs.nn) {
			if (v != 0) delete [] (v);
			nn=rhs.nn;
			v= new T[nn];
		}
		for (int i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}

template <class T>
Vec<T> & Vec<T>::operator=(const T &a)	//assign a to every element
{
	for (int i=0; i<nn; i++)
		v[i]=a;
	return *this;
}

template <class T>
inline T & Vec<T>::operator[](const int i)	//subscripting
{
	return v[i];
}

template <class T>
inline const T & Vec<T>::operator[](const int i) const	//subscripting
{
	return v[i];
}

template <class T>
inline int Vec<T>::size() const
{
	return nn;
}

template <class T>
T Vec<T>::min() const {
	if (nn == 0) {
		cout << "No element in the vector!!" << endl;
		return 0;
	} else if (nn == 1) {
		return v[0];
	} else {
		T min = v[0];
		for (int i=1; i<nn; i++) {
			if (v[i] < min) {
				min = v[i];
			}
		}
		return min;
	}
}

template <class T>
T Vec<T>::max() const {
	if (nn == 0) {
		cout << "No element in the vector!!" << endl;
		return 0;
	} else if (nn == 1) {
		return v[0];
	} else {
		T max = v[0];
		for (int i=1; i<nn; i++) {
			if (v[i] > max) {
				max = v[i];
			}
		}
		return max;
	}
}

template <class T>
void Vec<T>::print() const {
	if (nn == 0) {
		cout << "No element in the vector!!" << endl;
	} else {
		for (int i=0; i<nn; i++) {
			cout << v[i] << '\t';
		}
		cout << '\n';
	}	
}

template <class T>
void Vec<T>::write(ofstream &f) const {
	for (int i=0; i<nn; i++) {
		f << v[i] << '\n';
	}
}

template <class T>
void Vec<T>::write(string fname) const {
	ofstream out_file(fname.c_str());
	this->write(out_file);
	out_file.close();
}

template <class T>
T Vec<T>::dot_product(const Vec<T> &rhs) const {
	if (nn == 0) {
		cout << "No element in the vector!!" << endl;
		return 0;
	} else if (nn != rhs.nn) {
		cout << "The vector lengths for dot product are not equal!";
		return 0;
	} else {
		T prod = 0;
		for (int i=0; i<nn; i++) {
			prod += v[i]*rhs[i];
		}
		return prod;
	}
}

template <class T>
double Vec<T>::dist(const Vec<T> &rhs) const {
	double l = 0.0;
	if (nn == 0) {
		cout << "No element in the vector!!" << endl;
		return 0;
	} else if (nn != rhs.nn) {
		cout << "The vector lengths for dot product are not equal!";
		return 0;
	} else {
		for (int i=0; i<nn; i++) {
			l += (v[i]-rhs[i])*(v[i]-rhs[i]);
		}
		l = sqrt(l);
		return l;
	}
}

template <class T>
void Vec<T>::from_mat(const Mat<T> &mat, int dim, int id)
{
	if (dim == 1) {
		if (nn != mat.ncols()) {
			cout << "Length not equal!" << endl;
			return;
		}	
		for (int i=0; i<nn; i++) {
			v[i] = mat[id][i];
		}
	} else if (dim == 2) {
		if (nn != mat.nrows()) {
			cout << "Length not equal!" << endl;
			return;
		}	
		for (int i=0; i<nn; i++) {
			v[i] = mat[i][id];
		}
	} else {
		cout << "Error: parameter dim should be either 1:row or 2:col!" << endl;
	}
}

template <class T>
double Vec<T>::mean() const {
	double sum = 0.0;
	for (int i=0; i<nn; i++) {
		sum += v[i];
	}
	return (sum/nn);
}

template <class T>
double Vec<T>::var() const {
	double sum = 0.0;
	for (int i=0; i<nn; i++) {
		sum += v[i];
	}
	double mean = sum/nn;
	sum = 0.0;
	for (int i=0; i<nn; i++) {
		sum += (v[i]-mean)*(v[i]-mean);
	}
	return (sum/nn);
}

template <class T>
Vec<T>::~Vec()
{
	if (v != 0)
		delete[] (v);
}

template <class T>
class Mat {
private:
	int nn;
	int mm;
	T **v;
public:
	Mat();
	Mat(int n, int m);			// Zero-based array
	Mat(const T &a, int n, int m);	//Initialize to constant
	Mat(const T *a, int n, int m);	// Initialize to array
	Mat(ifstream &f, int n ,int m);	// Initialize to file
	Mat(const Mat &rhs);		// Copy constructor
	Mat & operator=(const Mat &rhs);	//assignment
	Mat & operator=(const T &a);		//assign a to every element
	inline T* operator[](const int i);	//subscripting: pointer to row i
	inline const T* operator[](const int i) const;
	inline int nrows() const;
	inline int ncols() const;
	void print() const;
	void write(ofstream &f) const;
	void write(string fname) const;

	~Mat();
};

template <class T>
Mat<T>::Mat() : nn(0), mm(0), v(0) {}

template <class T>
Mat<T>::Mat(int n, int m) : nn(n), mm(m), v(new T*[n])
{
	v[0] = new T[m*n];
	for (int i=1; i< n; i++)
		v[i] = v[i-1] + m;
}

template <class T>
Mat<T>::Mat(const T &a, int n, int m) : nn(n), mm(m), v(new T*[n])
{
	int i,j;
	v[0] = new T[m*n];
	for (i=1; i< n; i++)
		v[i] = v[i-1] + m;
	for (i=0; i< n; i++)
		for (j=0; j<m; j++)
			v[i][j] = a;
}

template <class T>
Mat<T>::Mat(const T *a, int n, int m) : nn(n), mm(m), v(new T*[n])
{
	int i,j;
	v[0] = new T[m*n];
	for (i=1; i< n; i++)
		v[i] = v[i-1] + m;
	for (i=0; i< n; i++)
		for (j=0; j<m; j++)
			v[i][j] = *a++;
}

template <class T>
Mat<T>::Mat(ifstream &f, int n, int m) : nn(n), mm(m), v(new T*[n])
{
	if (!f) {
		cout << "Error opening for construct Mat!!" << endl;
		return;
	}
	v[0] = new T[m*n];
	for(int i=1; i<n; i++)
		v[i] = v[i-1] + m;
	for(int i=0; i<n; i++)
		for( int j=0; j<m; j++)
			f >> v[i][j];
}


template <class T>
Mat<T>::Mat(const Mat &rhs) : nn(rhs.nn), mm(rhs.mm), v(new T*[nn])
{
	int i,j;
	v[0] = new T[mm*nn];
	for (i=1; i< nn; i++)
		v[i] = v[i-1] + mm;
	for (i=0; i< nn; i++)
		for (j=0; j<mm; j++)
			v[i][j] = rhs[i][j];
}

template <class T>
Mat<T> & Mat<T>::operator=(const Mat<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		int i,j;
		if (nn != rhs.nn || mm != rhs.mm) {
			if (v != 0) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.nn;
			mm=rhs.mm;
			v = new T*[nn];
			v[0] = new T[mm*nn];
		}
		for (i=1; i< nn; i++)
			v[i] = v[i-1] + mm;
		for (i=0; i< nn; i++)
			for (j=0; j<mm; j++)
				v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
Mat<T> & Mat<T>::operator=(const T &a)	//assign a to every element
{
	for (int i=0; i< nn; i++)
		for (int j=0; j<mm; j++)
			v[i][j] = a;
	return *this;
}

template <class T>
inline T* Mat<T>::operator[](const int i)	//subscripting: pointer to row i
{
	return v[i];
}

template <class T>
inline const T* Mat<T>::operator[](const int i) const
{
	return v[i];
}

template <class T>
inline int Mat<T>::nrows() const
{
	return nn;
}

template <class T>
inline int Mat<T>::ncols() const
{
	return mm;
}

template <class T>
void Mat<T>::print() const {
	for (int i=0; i<nn; i++) {
		for (int j=0; j<mm; j++) {
			cout << v[i][j] << '\t';
		}
		cout << '\n';
	}
}
template <class T>
void Mat<T>::write(ofstream &f) const {
	for (int i=0; i<nn; i++) {
		for (int j=0; j<mm; j++) {
			f << v[i][j] << '\t';
		}
		f << '\n';
	}
}
template <class T>
void Mat<T>::write(string fname) const {
	ofstream out_file(fname.c_str());
	this->write(out_file);
	out_file.close();
}


template <class T>
Mat<T>::~Mat()
{
	if (v != 0) {
		delete[] (v[0]);
		delete[] (v);
	}
}
#endif

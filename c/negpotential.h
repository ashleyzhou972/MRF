#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

double negpotential(double * y, double * ystar, int size_y, int ** neighbor,
	       	double alpha, double eta, double tau2);
double H_i(double * y, int i,  double * ystar, int size_y, int ** neighbor, 
		double alpha, double eta, double tau2);
double H_ij(double * y, int i, int j,  double * ystar, int size_y ,
		int ** neighbor, double alpha, double eta, double tau2);
double vector_multiplication(int * array1, double * array2, int size);


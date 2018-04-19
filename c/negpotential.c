#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>


typedef double (* pdf)(double, double *);
double negpotential(double* y, double * ystar, int size_y, double alpha, 
		double eta, double tau2); 
double H_ij(double y_i, double  ystar_i, double y_j, double ystar_j, pdf g);
double H_i(double y_i, double ystar_i, pdf g);
pdf mrf_normal(double x, double * params);



double negpotential( double * y, double * ystar, int size_y,  double alpha, double eta, double tau2, pdf g){
	double summand_i = 0;
	double summand_ij = 0;
	double result;
	for (int i = 0; i < size_y; i++){
		summand_i += H_i(y[i], ystar[i], g);
		for (int j = 0; j < size_y; j++){
			summand_ij += H_ij(y[i],ystar[i],y[j],ystar[j],g);
		}
	}
	result = summand_i + summand_ij;
	return result;
}

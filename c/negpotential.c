#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>


typedef double (* pdf)(double, double *);
double negpotential(double* w, double alpha, double eta, double tau2,); //TODO
double Hij(double y_i, double  ystar_i, double y_j, double ystar_j, 
	       pdf g);
double Hi(double y_i, double ystar_i, pdf g);



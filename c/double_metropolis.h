/**
 * @file double_metropolis.h
 * The double metropolis algorithm 
 * @author Naihui Zhou {nzhou@iastate.edu}
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <negpotential.h>


typedef double (* pdf)(double, double *);
//pdf functions, fisrt  x, second is pointer to other parameters;
typedef double (*negp) (double *, double);
//first x, second parameter. only one parameter;
typedef void (*auxiliary)(int, double *, double * int **, double *);
//generate auxiliary variable y, given special parameter;


/**
 * Main functions
 **/
double dm_step1(double theta0, pdf target, double var);
double dm_step2(int size_x,double * x, double theta_new, double theta_current, 
		int ** neighbor, negp negpotential, 
		auxiliary auxiliary_y_gibbs_theta);
/**
 * Distributional helper functions
 **/

double r_random_walk_chain(double current_y, double var);
double d_random_walk_chain(double current_y, double proposed_y,double var);
double jump_probability(double current, double proposed, pdf target);
void auxiliary_y_gibbs(int size_x, double * x, double * y, int ** neighbor, 
		double alpha, double eta, double tau2);
/**
 * Array helper functions
 **/
double vector_multiplication(int * array1, double * array2, int size);



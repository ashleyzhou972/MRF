/**
 * @file double_metropolis.h
 * The double metropolis algorithm 
 * @author Naihui Zhou {nzhou@iastate.edu}
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include "Rmath.h"
#include "negpotential.h"
#include "regular_metropolis.h"


typedef double (* pdf)(double, double *);
//pdf functions, fisrt  x, second is pointer to other parameters;
typedef double (*negp) (double *, double, double *, int **);
//first is array y; second parameter is the theta of interest;
//third is pointer to other parameters;
//fourth is neighbor matrix
typedef void (*auxiliary)(double *, double *,double *, int **);
//first is auxiliary array (y)
//second is starting array (x)
//third is pointer to other parameters 
//fourth is neighbor matrix


/**
 * Main functions
 **/
double dm_step1(double theta0, pdf target, double var, double *target_params);
double dm_step2(int size_x,double * x, double theta_new, double theta_current, 
		int ** neighbor, negp negpotential, 
		auxiliary auxiliary_y_gibbs_theta, double *par_neg, 
		double *par_auxi);
/**
 * Distributional helper functions
 **/

void auxiliary_y_gibbs(int size_x, double * x, double * y, int ** neighbor, 
		double alpha, double eta, double tau2);
/**
 * Array helper functions
 **/
double vector_multiplication(int * array1, double * array2, int size);

/**
 * Parameter-specific functions
 **/

double prior_alpha(double alpha_par, double *other_par);
double prior_eta(double eta_par, double *other_par);
double prior_tau2(double tau2_par, double *other_par);
double negp_alpha_t(double *y, double alpha_par, double *other_par, int ** neighbor);
double negp_eta_t(double *y, double eta_par, double *other_par, int ** neighbor);
double negp_tau2_t(double *y, double tau2_par, double *other_par, int ** neighbor);
void auxi_alpha_t(double *y,double *x, double *other_par, int **neighbor);
void auxi_eta_t(double *y,double *x, double *other_par, int **neighbor);
void auxi_tau2_t(double *y,double *x, double *other_par, int **neighbor);

double dm_step2_t_alpha(int t, double **w_bycol, double *alpha, double *eta, double *tau2, int N, int T, double alpha_new, int **neighbor);
double dm_step2_t_eta(int t, double **w_bycol, double *alpha, double *eta, double *tau2, int N, int T, double eta_new, int **neighbor);
double dm_step2_t_tau2(int t, double **w_bycol, double *alpha, double *eta, double *tau2, int N, int T, double tau2_new, int **neighbor);


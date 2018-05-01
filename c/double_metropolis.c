/**
 * 
 * Initial attempt at the double metropolis algorithm 
 * @author Naihui Zhou {nzhou@iastate.edu}
 * to run 
 * gcc -Wall -pedantic -o dm double_metropolis.c negpotential.c 
 * 	regular_metropolis.c -lm -lRmath 
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "negpotential.h"
#include "regular_metropolis.h"


typedef double (* pdf)(double, double *);
//pdf functions, fisrt  x, second is pointer to other parameters;
typedef double (*negp) (double *, double, double *, int **);
//first is array y; second parameter is the theta of interest;
//third is pointer to other parameters;
//fourth is neighbor matrix
typedef void (*auxiliary)(double *, double *, int **);
//first is auxiliary array
//second is pointer to other parameters 
//third is neighbor matrix


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

void auxiliary_y_gibbs(int size_x, double * x, double * y, int ** neighbor, 
		double alpha, double eta, double tau2);
/**
 * Array helper functions
 **/
double vector_multiplication(int * array1, double * array2, int size);





double dm_step1(double theta0, pdf target, double var, double *target_params){
	//in double-metropolis step 1,;
	//target_pdf is the prior distribution of theta;
	//target_pdf should return log scale
	double theta_new,jp,alpha;
	double u;
	do{
		theta_new = random_walk_chain(theta0, var);
		jp = jump_probability(theta0, theta_new, target, target_params);
		if (jp<1) alpha=jp;
		else alpha=1;
		u = runif(0.0,1.0);
	} while(u>alpha);
       return theta_new;	
}

double dm_step2(int size_x,double * x, double theta_new, double theta_current, 
		int ** neighbor, negp negpotential_theta, 
		auxiliary auxiliary_y_gibbs_theta){
	//target distribution here is the original 
	//data distribution f(x|theta) in the Bayes rule,
	//without the intractable constant;
	//aka the negpotential function;
	//The negpotential function should return in log scale;
	//First generate a auxiliary variable y;
	double * y[size_x];
	auxiliary_y_gibbs_theta(y);
	double r, alpha;
	double numerator, denominator;
	numerator = negpotential_theta(y,theta_current)+negpotential_theta(x, theta_new);
	denominator = negpotential_theta(x, theta_current)+negpotential_theta(y,theta_new);
	r = exp(numerator-denominator);
	if (r<1) alpha =r;
	else alpha = 1;
	double u;
	u = runif(0.0,1.0);
	if (u<=alpha) return theta_new;//accepted;
	else return theta_current;
}

void auxiliary_y_gibbs(int size_x, double * x, double * y, int ** neighbor, 
		double alpha, double eta, double tau2){
	//One Gibbs round for each element in x;
	//The generating pdf is the normal with mrf mu;
	int i =0 ;
	double mu_i, y_new_i;
	y = x;//starting value is x;
	//this may be redundant since y is initialized outside of this function;
	for (i =0;i<size_x;i++){
		mu_i = alpha + eta*vector_multiplication(neighbor[i], y,size_x);
		y_new_i = rnorm(mu_i, tau2);
		y[i] = y_new_i;
	}
}

double vector_multiplication(int * array1, double * array2, int size){
	double summand = 0;
	for (int i=0;i<size;i++){
		summand += (double) array1[i]*array2[i];
	}
	return summand;
}
/**
 * *******************************************************************
 * Above are the general functions
 * Below are wrappers for specific parameters, iterations, etc
 * ******************************************************************
 **/

/**
 * For alpha:
 **/
double prior_alpha(double alpha_par, double *other_par){
	//other_par should be the upper and lower bound of uniform distribution
	return dunif(alpha_par,other_par[0],other_par[1],1);
}

double negp_alpha_t(double *y, double alpha_par, double *other_par, int ** neighbor){
	//other_par should have N as first param
	//eta[t] as second param
	//tau2[t] as third param
	double ystar[N];
	for (int i=0;i<N;i++){
		ystar[i] = 0.0;
	}
	return negpotential(y,ystar,(int) other_par[0],neighbor,alpha_par,param[1],param[2]);
}


double auxi_alpha_t(double *y, double *x,double *other_par, int **neighbor){
	return auxiliary_y_gibbs((int)other_par[0], x, y, neighbor, other_par[1], other_par[2], other_par[3]);
}

double dm_step2_t_alpha(int t, double **w_bycol, double *alpha, double *eta, double *tau2, int N, int T, double alpha_new, int **neighbor){
	double params_alpha_prior[2];
	double params_alpha_neg[3];
	double params_alpha_aux[4];
	params_alpha_prior[0] = 0;
	params_alpha_neg[0] = N;
	
	return dm_step2(N, w_bycol[t], alpha_new, alpha[t], neighbor, negp_alpha_t, auxi_alpha_t);
}

/**
 * For eta:
 **/

double prior_eta(double eta_par, double *other_par){
	return dunif(eta-par, other_par[0], other_par[1], 1);
}

double negp_alpha_t(double *y, double alpha_par){







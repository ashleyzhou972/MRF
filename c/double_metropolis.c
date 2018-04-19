/**
 * 
 * Initial attempt at the double metropolis algorithm 
 * @author Naihui Zhou {nzhou@iastate.edu}
 * to run gcc -Wall -pedantic -o dm double_metropolis.c -lm -lRmath
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

typedef double (* pdf)(double, double *);
//pdf functions, fisrt  x, second is pointer to other parameters;
typedef double (*negp) (double *, double);
//first x, second parameter. only one parameter;


double data_pdf();
double mixing_pdf();
double prior_pdf();

double posterior_pdf_w();
/**
 * Each of the parameters alpha, eta and tau2
 * has two types of posteriors:
 * 	- The multivariate normal -1
 * 	- The negpotential function -2
 * For the negpotential function type posterior,
 * there are two methods to compute the intractable constant
 * 	- The double- metropolis algorithm
 * 	- importance sampling -3
 * Only case 2 will be in this script
 **/
double posterior_pdf_alpha();
double posterior_pdf_eta();
double posterior_pdf_tau2();

double dm_step1(double theta0, pdf target_pdf, double var);
double dm_step2(int size_x,double * x,pdf target, double theta_new, double theta_current);
/**
 * Assuming random walk chain means:
 *  	- proposal distribution q is normal with mean 0
 *  	- q(current|proposed) == q(proposed|current)
 *  	- all q terms get cancelled in calculating jump probability
 **/

double r_random_walk_chain(double current_y, double var){
	//simulate (and return) a new value from the proposal 
	//distribution (q), given current value;
	double z, y;
	z = rnorm(0.0,var);
	y = current_y + z;
	return y;
}


/**
 * this function is not used under random walk chain
 **/
double d_random_walk_chain(double current_y, double proposed_y,double var){
	//calculate density value of the proposal distribution (q)
	//should remain the same if exchange the current_y and proposed_y
	//positions
	return dnorm(current_y-proposed_y,0.0,var);
}

double jump_probability(double current, double proposed, pdf target){
	//assuming random walk chain
	double alpha;
	double numerator, denominator;
	double * dummy_param[2];
	numerator = target(proposed, dummy_param); 
	denominator = target(current, dummy_param); 
	alpha = numerator/denominator;
	return alpha;
}
double dm_step1(double theta0, pdf target_pdf, double var){
	//in double-metropolis step 1,;
	//target_pdf is the prior distribution of theta;
	double theta_new,jp,alpha;
	double u;
	do{
		theta_new = random_walk_chain(theta0, var);
		jp = jump_probability(theta0, theta_new, target_pdf);
		if (jp<1) alpha=jp;
		else alpha=1;
		u = runif(0.0,1.0);
	} while(u>alpha);
       return theta_new;	
}

double dm_step2(int size_x,double * x,pdf target, double theta_new, double theta_current){
	//target distribution here is the original 
	//data distribution f(x|theta) in the Bayes rule,
	//without the intractable constant;
	//aka the negpotential function;
	//First generate a auxiliary variable y;
	double * y[size_x];
	//TODO Gibbs procedure to simulate the auxiliary y;
	double r, alpha;
	double numerator, denominator;
	numerator = negpotential(y,theta_current)*negpotential(x, theta_new);
	denominator = negpotential(x, theta_current)*negpotential(y,theta_new);
	r = numerator/denominator;
	if (r<1) alpha =r;
	else alpha = 1;
	double u;
	u = runif(0.0,1.0);
	if (u<=alpha) return theta_new;//accepted;
	else return theta_current;
}

void auxiliary_y_gibbs(double * x, double * y, _Bool ** neighbor, double tau2){

}


int main(void){
	return 0;
}

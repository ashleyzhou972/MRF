#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "regular_metropolis.h"
#include "negpotential.h"

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

/**
 * Assuming random walk chain means:
 *  	- proposal distribution q is normal with mean 0
 *  	- q(current|proposed) == q(proposed|current)
 *  	- all q terms get cancelled in calculating jump probability
 **/

double r_random_walk_chain(double current_y, double var)
{
	//simulate (and return) a new value from the proposal
	//distribution (q), given current value;
	double z, y;

	z = rnorm(0.0, var);
	y = current_y + z;

	return y;
}


/**
 * ************************************
 * this function is not used under random walk chain
 **/
double d_random_walk_chain(double current_y, double proposed_y, double var)
{
	//calculate density value of the proposal distribution (q)
	//should remain the same if exchange the current_y and proposed_y
	//positions
	return dnorm(current_y-proposed_y, 0.0, var, 1);
}
/**
 **************************************
 **/

/**
 * Assuming random walk chain means:
 *  	- proposal distribution q is normal with mean 0
 *  	- q(current|proposed) == q(proposed|current)
 *  	- all q terms get cancelled in calculating jump probability
 **/


double jump_probability(double current, double proposed, pdf target,
		double *param)
{
	//assuming random walk chain
	//Output of target should be log scaled
	double alpha, prob;
	double numerator, denominator;

	numerator = target(proposed, param);
	denominator = target(current, param);
	alpha = numerator-denominator;
	prob = exp(alpha);
	if (prob < 1.0)
		return prob;
	else
		return 1.0;
}


/**
 * The log density of univariate y given w
 * This function drops the factorial!!
 * Because factorials overflow
 * (1/y!) here can be seen as a constant with regard to w in the posterior
 **/
double log_data_density_univar(double y, double w)
{
	return -exp(w) + w*y;
}

//w is the column observation at w[t];	
//w_in may not be the same as w[i];
//w_in can be newly generated in the Metropolis process;
double log_mrf_density_univar(int size_w, int i, double w_in, double *w, int **neighbor, double alpha, double eta, double tau2)
{
	double mu_i;
	double w_subtracted[size_w];
	for (int k = 0; k<size_w; ++k) 
		w_subtracted[k] = w[k] - alpha;
	mu_i = alpha + eta*vector_multiplication(neighbor[i], w_subtracted, size_w);
	return dnorm(w_in, mu_i, tau2, 1);
}

double log_sum_density_univar(int size_w, int i, double y, double w_in, double *w, int **neighbor, double alpha, double eta, double tau2)
{
	double sum;
	sum = log_data_density_univar(y, w_in) + log_mrf_density_univar(size_w, i, w_in, w, neighbor, alpha, eta, tau2);
	return sum;
}


/**
 * @param t  the current iteration
 * try to fill (t+1)th iteration
 * @param w N by T matrix as a double pointer
 * @param y observed data array
 * @param var variance for simulating new values in metropolis
 * @param data_pdf log pdf of data distribution f(y_i|w_i)
 **/
int metropolis_for_w_univar(int t, int N, double **w, double *y, double var, int **neighbor, double alpha, double eta, double tau2)
{
	// w is a N by T matrix;
	// y is a size N vector
	int jumps = 0;
	for (int i = 0; i < N; ++i) {
		double w_new_i, prob, jp;
		double numerator, denominator;
		double u = runif(0.0, 1.0);

		w_new_i = r_random_walk_chain(w[t][i], var);
		/**
		 * calculating jump probability
		 * not using the jump_probability function
		 * updated 20180523
		 **/
		numerator = log_sum_density_univar(N, i, y[i], w_new_i, w[t], neighbor, alpha, eta, tau2);
		denominator = log_sum_density_univar(N, i, y[i], w[t][i], w[t], neighbor, alpha, eta, tau2);

		prob = exp(numerator - denominator );
		if (prob < 1.0)
			jp = prob;
		else
			jp = 1;
		//end of calculating jump probability;
		if (u < jp) {
			w[t+1][i] = w_new_i;
			jumps += 1;
		}
		else 
			w[t+1][i] = w[t][i];
	}
	return jumps;
}

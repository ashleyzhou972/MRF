#include <stdio.h>
#include <stdlib.h>

typedef double (* pdf)(double, double *);
//pdf functions, fisrt  x, second is pointer to other parameters;
typedef double (*negp) (double *, double);
//first x, second parameter. only one parameter;
typedef void (*auxiliary)(int, double *, double *,int **, double *);
//generate auxiliary variable y, given special parameter;



/**
 *
 * For regular metropolis-hastings within Gibbs algorithm
 **/
double r_random_walk_chain(double current_y, double var);
double d_random_walk_chain(double current_y, double proposed_y,double var);
double jump_probability(double current, double proposed, pdf target);

/**
 * regular metropolis for w
 **/

double log_data_density_univar(double y, double w);



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

double r_random_walk_chain(double current_y, double var){
	//simulate (and return) a new value from the proposal 
	//distribution (q), given current value;
	double z, y;
	z = rnorm(0.0,var);
	y = current_y + z;
	return y;
}


/**
 * ************************************
 * this function is not used under random walk chain
 **/
double d_random_walk_chain(double current_y, double proposed_y,double var){
	//calculate density value of the proposal distribution (q)
	//should remain the same if exchange the current_y and proposed_y
	//positions
	return dnorm(current_y-proposed_y,0.0,var);
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
		double *param){
	//assuming random walk chain
	//Output of target should be log scaled
	double alpha, prob;
	double numerator, denominator;
	numerator = target(proposed, param); 
	denominator = target(current, param); 
	alpha = numerator-denominator;
	prob = exp(alpha);
	if (prob<1) return prob
	else return 1
}


/**
 * 
 * The log density of univariate y given w
 * This function drops the factorial!! 
 * Because factorials overflow 
 * (1/y!) here can be seen as a constant with regard to w in the posterior
 **/
double log_data_density_univar(double y, double w){
	return -exp(w) + w*y;
}

double pdf_lddu(double w, double *y){
	return log_data_density_univar(*y, w);
}

/**
 * @param t  the current iteration
 * try to fill (t+1)th iteration
 * @param w N by T matrix as a double pointer
 * @param y observed data array
 * @param var variance for simulating new values in metropolis
 * @param data_pdf log pdf of data distribution f(y_i|w_i)
 **/
void metropolis_for_w_univar(int t, double **w, double *y, double var){
	// w is a N by T matrix;
	// y is a size N vector
	for (i=0;i<N;i++){
		double w_new_i;
		w_new_i = r_random_walk_chain(w[i][t],var);
		jp = jump_probability(w[i][t],w_new_i,pdf_lddu, &y[i]);
		double u = runif(0,1);
		if (u<jp) w[i][t+1] = w_new;
		else w[i][t+1] = w[i][t];
	}
}

		


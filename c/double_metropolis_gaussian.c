/**
 * @file double_metropolis.c
 * Initial attempt at the double metropolis algorithm
 * @author Naihui Zhou {nzhou@iastate.edu}
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include "Rmath.h"
#include "negpotential.h"
#include "regular_metropolis.h"
#include "double_metropolis_gaussian.h"


/**
 * The Gaussian case is when we have OBSERVED gaussian spatial MRF
 * as opposed to an observed Poisson data and mixed with a Gaussian MRF
 **/

double dm_step1(double theta0, pdf target, double var, double *target_params)
{
	//in double-metropolis step 1,;
	//target_pdf is the prior distribution of theta;
	//target_pdf should return log scale
	double theta_new, jp, alpha;
	double u;

	theta_new = r_random_walk_chain(theta0, var);
	jp = jump_probability(theta0, theta_new, target, target_params);
	if (jp < 1.0)
		alpha = jp;
	else
		alpha = 1.0;
	u = runif(0.0, 1.0);
	if (u < alpha)
		return theta_new;
	else
		return theta0;
}

int dm_step2(int size_x, double *x, double theta_new, double theta_current,
		int **neighbor, negp negpotential_theta,
		auxiliary auxiliary_y_gibbs_theta, double *par_neg,
		double *par_auxi)
{
	//target distribution here is the original
	//data distribution f(x|theta) in the Bayes rule,
	//without the intractable constant;
	//aka the negpotential function;
	//The negpotential function should return in log scale;
	//First generate a auxiliary variable y;
	double y[size_x];
	double r, alpha, u;
	double numerator, denominator;

	auxiliary_y_gibbs_theta(y, x, par_auxi, neighbor);
	numerator = negpotential_theta(y, theta_current, par_neg, neighbor)
		+ negpotential_theta(x, theta_new, par_neg, neighbor);

	denominator = negpotential_theta(x, theta_current, par_neg, neighbor)
		+ negpotential_theta(y, theta_new, par_neg, neighbor);

	r = exp(numerator-denominator);
	if (r < 1)
		alpha = r;
	else
		alpha = 1;

	u = runif(0.0, 1.0);
	if (u <= alpha)
		return 1;//accepted;
	else
		return 0;
}

void auxiliary_y_gibbs(int size_x, double *x, double *y, int **neighbor,
		double alpha, double eta, double tau2)
{
	//One Gibbs round for each element in x;
	//The generating pdf is the normal with mrf mu;
	int i = 0;
	double mu_i, y_new_i;

	for (int k = 0; k < size_x; ++k)
		y[k] = x[k]; //starting value is x

	//this may be redundant since y is initialized outside of this function;
	for (i = 0; i < size_x; ++i) {
		mu_i = alpha + eta*vector_multiplication(neighbor[i], y, size_x);
		y_new_i = rnorm(mu_i, tau2);
		y[i] = y_new_i;
	}
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
double prior_alpha(double alpha_par, double *other_par)
{
	//other_par should be the upper and lower bound of uniform distribution
	double lb = other_par[0];
	double ub = other_par[1];

	return dunif(alpha_par, lb, ub, 1);
}

double negp_alpha_t(double *y, double alpha_par, double *other_par, int **neighbor)
{
	//other_par should have N as first param
	//eta[t] as second param
	//tau2[t] as third param
	int N = (int) other_par[0];
	double eta = other_par[1];
	double tau2 = other_par[2];
	double ystar[N];

	for (int i = 0; i < N; ++i)
		ystar[i] = 0.0;

	return negpotential(y, ystar, N, neighbor, alpha_par, eta, tau2);
}


void auxi_alpha_t(double *y, double *x, double *other_par, int **neighbor)
{
	int N = (int) other_par[0];
	double alpha_new = other_par[1];//!this is the new alpha generated;
	double eta = other_par[2];
	double tau2 = other_par[3];

	auxiliary_y_gibbs(N, x, y, neighbor, alpha_new, eta, tau2);
}

int dm_step2_t_alpha(int t, double *y, double *alpha, double *eta,
		double *tau2, int N, int T, double alpha_new, int **neighbor)
{
	double par_neg[3];
	double par_auxi[4];

	par_neg[0] = N;
	par_neg[1] = eta[t];
	par_neg[2] = tau2[t];
	par_auxi[0] = N;
	par_auxi[1] = alpha_new;
	par_auxi[2] = eta[t];
	par_auxi[3] = tau2[t];

	return dm_step2(N, y, alpha_new, alpha[t], neighbor,
			negp_alpha_t, auxi_alpha_t, par_neg, par_auxi);
}

/**
 * For eta:
 **/
double prior_eta(double eta_par, double *other_par)
{
	//other_par should be the upper and lower bound of uniform distribution
	double lb = other_par[0];
	double ub = other_par[1];

	return dunif(eta_par, lb, ub, 1);
}

double negp_eta_t(double *y, double eta_par, double *other_par, int **neighbor)
{
	//other_par should have N as first param
	//alpha[t+1] as second param
	//tau2[t] as third param
	int N = (int) other_par[0];
	double ystar[N];
	double alpha = other_par[1];
	double tau2 = other_par[2];

	for (int i = 0; i < N; ++i)
		ystar[i] = 0.0;

	return negpotential(y, ystar, N, neighbor, alpha, eta_par, tau2);
}


void auxi_eta_t(double *y, double *x, double *other_par, int **neighbor)
{
	int N = (int) other_par[0];
	double alpha = other_par[1];
	double eta_new = other_par[2];//This is the new eta from dm_step1
	double tau2 = other_par[3];

	auxiliary_y_gibbs(N, x, y, neighbor, alpha, eta_new, tau2);

}

int dm_step2_t_eta(int t, double *y, double *alpha, double *eta,
		double *tau2, int N, int T, double eta_new, int **neighbor)
{
	double par_neg[3];
	double par_auxi[4];

	par_neg[0] = N;
	par_neg[1] = alpha[t+1];
	par_neg[2] = tau2[t];
	par_auxi[0] = N;
	par_auxi[1] = alpha[t+1];
	par_auxi[2] = eta_new;
	par_auxi[3] = tau2[t];

	return dm_step2(N, y, eta_new, eta[t], neighbor, negp_eta_t,
			auxi_eta_t, par_neg, par_auxi);
}


/**
 * For tau2
 **/
double prior_tau2(double tau2_par, double *other_par)
{
	//other_par should be the upper and lower bound of uniform distribution
	double lb = other_par[0];
	double ub = other_par[1];

	return dunif(tau2_par, lb, ub, 1);
}

double negp_tau2_t(double *y, double tau2_par, double *other_par, int **neighbor)
{
	//other_par should have N as first param
	//alpha[t+1] as second param
	//eta[t+1] as third param
	int N = (int) other_par[0];
	double ystar[N];
	double alpha = other_par[1];
	double eta = other_par[2];

	for (int i = 0; i < N; ++i)
		ystar[i] = 0.0;

	return negpotential(y, ystar, N, neighbor, alpha, eta, tau2_par);
}


void auxi_tau2_t(double *y, double *x, double *other_par, int **neighbor)
{
	int N = (int) other_par[0];
	double alpha = other_par[1];
	double eta = other_par[2];
	double tau2_new = other_par[3];//!this is the new tau2 from dm_step1

	auxiliary_y_gibbs(N, x, y, neighbor, alpha, eta, tau2_new);
}

int dm_step2_t_tau2(int t, double *y, double *alpha, double *eta,
		double *tau2, int N, int T, double tau2_new, int **neighbor)
{
	double par_neg[3];
	double par_auxi[4];

	par_neg[0] = N;
	par_neg[1] = alpha[t+1];
	par_neg[2] = eta[t+1];
	par_auxi[0] = N;
	par_auxi[1] = alpha[t+1];
	par_auxi[2] = eta[t+1];
	par_auxi[3] = tau2_new;

	return dm_step2(N, y, tau2_new, tau2[t], neighbor, negp_tau2_t,
			auxi_tau2_t, par_neg, par_auxi);
}

/**
 * How to compile
 * gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG      -fpic  -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c dm_call.c -o dm_call.o -fopenmp -DOPENMP
 * gcc -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions -Wl,-z,relro -o dm_call.so dm_call.o double_metropolis.o regular_metropolis.o negpotential.o -lm -lRmath -L/usr/lib/R/lib -lR -lRmath
 **/

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*#define MATHLIB_STANDALONE*/
#include <Rmath.h>
#include "double_metropolis_gaussian.h"
#include "regular_metropolis.h"
#include "negpotential.h"

/**
 * The Gaussian case is when we have OBSERVED gaussian spatial MRF
 * as opposed to an observed Poisson data and mixed with a Gaussian MRF
 **/

void initialize(double alpha0, double eta0, double tau20,
		double *alpha, double *eta, double *tau2, int N, int T);

SEXP double_metropolis_gaussian(SEXP T_in, SEXP y_in, SEXP neighbor_in, SEXP vars_in, SEXP bounds_alpha, SEXP bounds_eta, SEXP bounds_tau2, SEXP initials)
{
	//declaration of variables;
	int N, T;
	int *T_p;
/*	int *y_int;*/
	double *y;
	int *neighbor_1d;
	int dim_neighbor;
	int num_protected = 0;
	double *vars;
	double *b_alpha, *b_eta, *b_tau2;
	double *initials_temp;
	double alpha0, eta0, tau20;
	//Check incoming SEXP types and extract data
	if (!isInteger(T_in))
		error("[ERROR] First argument must be interger");
	else {
		T_p = INTEGER(T_in);
		T = *T_p;
	}
	if (!isReal(y_in) || !isVector(y_in))
		error("[ERROR] Second argument must be integer vector");
	else {
		y = REAL(y_in);
		N = length(y_in);
	}

	int **neighbor_2d;

	neighbor_2d = malloc(N*sizeof(int *));
	for (int i = 0; i < N; ++i)
		neighbor_2d[i] = malloc(N*sizeof(int));
	if (!isInteger(neighbor_in) || !isVector(neighbor_in))
		error("[ERROR] Third argument must be integer vector");
	else {
		neighbor_1d = INTEGER(neighbor_in);
		dim_neighbor = length(neighbor_in);
		if (dim_neighbor != N*N)
			error("[ERROR] Size of neighborhood matrix does not agree with size of data vector");
		else {
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j)
					neighbor_2d[i][j] = neighbor_1d[j*N+i];
			}
		}
	}
	if (!isReal(vars_in) || !isVector(vars_in) || !length(vars_in) == 4)
		error("[ERROR] Fourth argument must be real vector with 4 elements");
	else
		vars = REAL(vars_in);
	if (!isReal(bounds_alpha) || !isReal(bounds_eta) || !isReal(bounds_tau2)
			|| !isVector(bounds_alpha) || !isVector(bounds_eta) ||
			!isVector(bounds_tau2))
		error("[ERROR] three bounds arguments must be real vectors\n");
	else {
		if (!length(bounds_alpha) == 2 || !length(bounds_eta) == 2 ||
					!length(bounds_tau2) == 2)
			error("[ERROR] three bounds arguments must each be of length 2\n");

		else {
			b_alpha = REAL(bounds_alpha);
			b_eta = REAL(bounds_eta);
			b_tau2 = REAL(bounds_tau2);
		}
	}
	if (!isReal(initials) || !isVector(initials))
		error("[ERROR] last argument should be a real vector\n");
	else {
		initials_temp = REAL(initials);
		alpha0 = initials_temp[0];
		eta0 = initials_temp[1];
		tau20 = initials_temp[2];
	}
	Rprintf("size of data vector is %d\n", N);
	Rprintf("number of iterations is %d\n", T);

	SEXP R_alpha, R_eta, R_tau2, R_Return_List;

	PROTECT(R_alpha = allocVector(REALSXP, T+1));
	PROTECT(R_eta = allocVector(REALSXP, T+1));
	PROTECT(R_tau2 = allocVector(REALSXP, T+1));
	PROTECT(R_Return_List = allocVector(VECSXP, 4));
	num_protected += 4;

/*	allocate(&w, N, T+1);*/
	initialize(alpha0, eta0, tau20, REAL(R_alpha), REAL(R_eta), REAL(R_tau2), N, T+1);

	int t; //iteration counter;

	/**
	 *	within iteration t:
	 *	- step1 Use dm (double metropolis) to generate posterior alpha
	 *	- step2 Use dm to generate posterior eta
	 *	- step3 Use dm to generate posterior tau2
	 **/
	double new_alpha, new_eta, new_tau2;
	double *alpha = REAL(R_alpha);
	double *eta = REAL(R_eta);
	double *tau2 = REAL(R_tau2);

	int jc_alpha = 0;
	int jc_eta = 0;
	int jc_tau2 = 0;
	int ret_alpha, ret_eta, ret_tau2;
	for (t = 0; t < T; ++t) {
		//printf("MC Iteration %d\n", t+1);
		//step1 (alpha);
		new_alpha = dm_step1(alpha[t], prior_alpha, vars[1], b_alpha);
		ret_alpha = dm_step2_t_alpha(t, y, alpha, eta, tau2, N, T, new_alpha, neighbor_2d);
		if (ret_alpha == 1) {
			alpha[t+1] = new_alpha;
			jc_alpha += 1;
		}
		else 
			alpha[t+1] = alpha[t];
		//step2 (eta);
		new_eta = dm_step1(eta[t], prior_eta, vars[2], b_eta);
		ret_eta = dm_step2_t_eta(t, y, alpha, eta, tau2, N, T, new_eta, neighbor_2d);
		if (ret_eta == 1) {
			eta[t+1] = new_eta;
			jc_eta += 1;
		}
		else 
			eta[t+1] = eta[t];
		//step3 (tau2);
		new_tau2 = dm_step1(tau2[t], prior_tau2, vars[3], b_tau2);
		ret_tau2 = dm_step2_t_tau2(t, y, alpha, eta, tau2, N, T, new_tau2, neighbor_2d);
		if (ret_tau2 == 1) {
			tau2[t+1] = new_tau2;
			jc_tau2 += 1;
		}
		else
			tau2[t+1] = tau2[t];
	}
	//printf("jump counts are %d, %d, %d, %d\n", jc_w, jc_alpha, jc_eta, jc_tau2);
	SET_VECTOR_ELT(R_Return_List, 0, R_alpha);
	SET_VECTOR_ELT(R_Return_List, 1, R_eta);
	SET_VECTOR_ELT(R_Return_List, 2, R_tau2);
	//begin ;
	SEXP R_jc;
	PROTECT(R_jc = allocVector(INTSXP,3));
	num_protected += 1;
	int *jc = INTEGER(R_jc);
	jc[0] = jc_alpha;
	jc[1] = jc_eta;
	jc[2] = jc_tau2;
	SET_VECTOR_ELT(R_Return_List, 3, R_jc);
	//end;
	
	UNPROTECT(num_protected);
	return R_Return_List;
}



void initialize(double alpha0, double eta0, double tau20,
		double *alpha, double *eta, double *tau2, int N, int T)
{
	alpha[0] = alpha0;
	eta[0] = eta0;
	tau2[0] = tau20;
}

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
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "cblas_double_metropolis.h"
#include "ms_regular_metropolis.h"
#include "cblas_negpotential.h"

void allocate_column(double *w, double **w_bycol, int N, int T);
void initialize_from_input(double *w_in, double alpha0, double eta0, 
		double tau20, double *w, double *alpha, double *eta, 
		double *tau2, int N, int T);
SEXP double_metropolis_cont(SEXP T_in, SEXP N_in, SEXP y_in, SEXP neighbor_in, SEXP vars_in, SEXP bounds_alpha, SEXP bounds_eta, SEXP bounds_tau2, SEXP initials, SEXP wInitials)
{
	//declaration of variables;
	//m is number of samples;
	int N, T, m;
	int *T_p;
	int *N_p;
/*	int *y_int;*/
	double *y;
	int *neighbor_1d;
	int dim_neighbor;
	int num_protected = 0;
	double *vars;
	double *b_alpha, *b_eta, *b_tau2;
	double *initials_temp;
	double *w_ini;
	double alpha0, eta0, tau20;
	//Check incoming SEXP types and extract data
	if (!isInteger(T_in))
		error("[ERROR] Iteration argument must be integer");
	else {
		T_p = INTEGER(T_in);
		T = *T_p;
	}
	if (!isInteger(N_in))
		error("[ERROR] Size argument must be integer");
	else {
		N_p = INTEGER(N_in);
		N = *N_p;
	}
	if (!isReal(y_in) || !isVector(y_in))
		error("[ERROR] y argument must be real vector");
	else {
		y = REAL(y_in); //This is m samples collapsed by column;
		m = length(y_in) / N;
		printf("Input %d samples\n", m);
	}

	int **neighbor_2d;

	neighbor_2d = malloc(N*sizeof(int *));
	for (int i = 0; i < N; ++i)
		neighbor_2d[i] = malloc(N*sizeof(int));
	if (!isInteger(neighbor_in) || !isVector(neighbor_in))
		error("[ERROR] Neighborhood argument must be integer vector");
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
/**
 * updated 20180828
 * take initial w from R function input
 * can be the last w from previous MCMC
**/
	if (!isReal(wInitials) || !isVector(wInitials))
		error("[ERROR] initial w argument should be a real vector\n");
	else {
		w_ini = REAL(wInitials);
	}

	Rprintf("size of data vector is %d\n", N);
	Rprintf("number of iterations is %d\n", T);

	SEXP R_alpha, R_eta, R_tau2, R_w, R_Return_List;

	PROTECT(R_w = allocMatrix(REALSXP, N, T+1));
	PROTECT(R_alpha = allocVector(REALSXP, T+1));
	PROTECT(R_eta = allocVector(REALSXP, T+1));
	PROTECT(R_tau2 = allocVector(REALSXP, T+1));
	PROTECT(R_Return_List = allocVector(VECSXP, 5));
	num_protected += 5;

/*	allocate(&w, N, T+1);*/
	initialize_from_input(w_ini, alpha0, eta0, tau20, REAL(R_w), REAL(R_alpha), REAL(R_eta), REAL(R_tau2), N, T+1);

	int t; //iteration counter;

	/**
	 *	within iteration t:
	 *	- step1 Use rm (regular metropolis) to generate posterior w;
	 *	- step2 Use dm (double metropolis) to generate posterior alpha
	 *	- step2 Use dm to generate posterior eta
	 *	- step3 Use dm to generate posterior tau2
	 **/
	double new_alpha, new_eta, new_tau2;
	double *alpha = REAL(R_alpha);
	double *eta = REAL(R_eta);
	double *tau2 = REAL(R_tau2);
	double *w = REAL(R_w);
	double **w_bycol = malloc(T*sizeof **w_bycol);

	int jc_w = 0;
	int jc_alpha = 0;
	int jc_eta = 0;
	int jc_tau2 = 0;
	int ret_w, ret_alpha, ret_eta, ret_tau2;
	allocate_column(w, w_bycol, N, T+1);
	for (t = 0; t < T; ++t) {
		if (t%10==0) {
			fprintf(stdout, "MC Iteration %d\n", t+1);
			//printf("MC Iteration %d\n", t+1);
		}		
		//step1;
		ret_w = metropolis_for_w_vector_mu(t, N, w_bycol, y, m,  vars[0], neighbor_1d, alpha[t], eta[t], tau2[t]);
		jc_w += ret_w;
		//step2 (alpha);
		new_alpha = dm_step1(alpha[t], prior_alpha, vars[1], b_alpha);
		ret_alpha = dm_step2_t_alpha(t, w_bycol, alpha, eta, tau2, N, T, new_alpha, neighbor_2d);
		if (ret_alpha == 1) { //accepted;
			alpha[t+1] = new_alpha;
			jc_alpha += 1;
		}
		else 
			alpha[t+1] = alpha[t];
		//step3 (eta);
		new_eta = dm_step1(eta[t], prior_eta, vars[2], b_eta);
		ret_eta = dm_step2_t_eta(t, w_bycol, alpha, eta, tau2, N, T, new_eta, neighbor_2d);
		if (ret_eta == 1) {
			eta[t+1] = new_eta;
			jc_eta += 1;
		}
		else 
			eta[t+1] = eta[t];
		//step4 (tau2);
		new_tau2 = dm_step1(tau2[t], prior_tau2, vars[3], b_tau2);
		ret_tau2 = dm_step2_t_tau2(t, w_bycol, alpha, eta, tau2, N, T, new_tau2, neighbor_2d);
		if (ret_tau2 == 1) {
			tau2[t+1] = new_tau2;
			jc_tau2 += 1;
		}
		else
			tau2[t+1] = tau2[t];
	}
	printf("jump counts are %d, %d, %d, %d\n", jc_w, jc_alpha, jc_eta, jc_tau2);
	SET_VECTOR_ELT(R_Return_List, 0, R_w);
	SET_VECTOR_ELT(R_Return_List, 1, R_alpha);
	SET_VECTOR_ELT(R_Return_List, 2, R_eta);
	SET_VECTOR_ELT(R_Return_List, 3, R_tau2);
	//begin ;
	SEXP R_jc;
	PROTECT(R_jc = allocVector(INTSXP,4));
	num_protected += 1;
	int *jc = INTEGER(R_jc);
	jc[0] = jc_w;
	jc[1] = jc_alpha;
	jc[2] = jc_eta;
	jc[3] = jc_tau2;
	SET_VECTOR_ELT(R_Return_List, 4, R_jc);
	//end;
	
	UNPROTECT(num_protected);
	return R_Return_List;
}


void allocate_column(double *w, double **w_bycol, int N, int T)
{
	for (int i = 0; i < T; ++i)
		w_bycol[i] = (double *)w + i*N;
}

/**
 * w is initialized to have rnorm(2,1)
 * w has N rows and T columns
 * w[i,j] = w[i*N+j]
 * we fill out the first column (iteration 0)
 **/

void initialize_from_input(double *w_in, double alpha0, double eta0, 
		double tau20, double *w, double *alpha, double *eta, 
		double *tau2, int N, int T)
{
	alpha[0] = alpha0;
	eta[0] = eta0;
	tau2[0] = tau20;
	for (int i =0; i< N; i++) {
		w[i] = w_in[i];
	}
}

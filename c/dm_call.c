#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*#define MATHLIB_STANDALONE*/
#include <Rmath.h>
#include "double_metropolis.h"
#include "regular_metropolis.h"
#include "negpotential.h"

void allocate_column(double *w, double **w_bycol, int N, int T);
void initialize(double alpha0, double eta0, double tau20, double *w,
		double *alpha, double *eta, double *tau2, int N, int T);

SEXP double_metropolis(SEXP T_in, SEXP y_in, SEXP neighbor_in, SEXP vars_in, SEXP bounds_alpha, SEXP bounds_eta, SEXP bounds_tau2, SEXP initials){
	//declaration of variables;
	int N,T;
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
		if (dim_neighbor!=N*N)
			error("[ERROR] Size of neighborhood matrix does not agree with size of data vector");
		else {
			for (int i = 0;i<N;i++) {
				for (int j = 0;j<N;j++) 
					neighbor_2d[i][j] = neighbor_1d[j*N+i];
			}
		}
	}
	if (!isReal(vars_in) || !isVector(vars_in) || !length(vars_in)==4)
		error("[ERROR] Fourth argument must be real vector with 4 elements");
	else 
		vars = REAL(vars_in);
	if (!isReal(bounds_alpha) || !isReal(bounds_eta) || !isReal(bounds_tau2) 
			|| !isVector(bounds_alpha) || !isVector(bounds_eta) ||
			!isVector(bounds_tau2)) 
		error("[ERROR] three bounds arguments must be real vectors\n");
	else {
		if (!length(bounds_alpha)==2 || !length(bounds_eta)==2 ||
					!length(bounds_tau2)==2)
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
	Rprintf("size of data vector is %d\n" , N);
	Rprintf("number of iterations is %d\n" , T);
	/*
	Rprintf("data vector is :\n");
	for (int i =0;i<N;i++) {
		Rprintf("%d " , y_int[i]);
	}
	Rprintf("\n");
	Rprintf("first column of neighbor matrix is :\n");
	for (int i=0;i<N;i++){
		Rprintf("%d ",neighbor_2d[i][0]);
	}
	Rprintf("\n");
	*/
	
	SEXP R_alpha, R_eta, R_tau2, R_w, R_Return_List;
	PROTECT(R_w = allocMatrix(REALSXP, N, T+1));
	PROTECT(R_alpha = allocVector(REALSXP, T+1));
	PROTECT(R_eta = allocVector(REALSXP, T+1));
	PROTECT(R_tau2 = allocVector(REALSXP, T+1));
	PROTECT(R_Return_List = allocVector(VECSXP, 3));
	num_protected += 5;
	

/*	allocate(&w, N, T+1);*/
	initialize(alpha0, eta0, tau20, REAL(R_w), REAL(R_alpha), REAL(R_eta), REAL(R_tau2), N, T+1);

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
	allocate_column(w, w_bycol, N, T+1);
	for (t = 0; t < T; ++t) {
		printf("MC Iteration %d\n", t+1);
		//step1;
		metropolis_for_w_univar(t, N, w_bycol, y, vars[0]);
		//step2 (alpha);
		new_alpha = dm_step1(alpha[t], prior_alpha, vars[1],b_alpha);
		alpha[t+1] = dm_step2_t_alpha(t, w_bycol, alpha, eta, tau2, N, T, new_alpha, neighbor_2d);
		//step3 (eta);
		new_eta = dm_step1(eta[t], prior_eta, vars[2], b_eta);
		eta[t+1] = dm_step2_t_eta(t, w_bycol, alpha, eta, tau2, N, T, new_eta, neighbor_2d);
		//step4 (tau2);
		new_tau2 = dm_step1(tau2[t], prior_tau2, vars[3], b_tau2);
		tau2[t+1] = dm_step2_t_tau2(t, w_bycol, alpha, eta, tau2, N, T, new_tau2, neighbor_2d);
	}
	free(w_bycol[0]);
	free(w_bycol);
	for (int k = 0; k < N; ++k)
		free(neighbor_2d[k]);
	free(neighbor_2d);

	SET_VECTOR_ELT(R_Return_List, 0, R_w);
	SET_VECTOR_ELT(R_Return_List, 1, R_alpha);
	SET_VECTOR_ELT(R_Return_List, 2, R_eta);
	SET_VECTOR_ELT(R_Return_List, 3, R_tau2);
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
void initialize(double alpha0, double eta0, double tau20, double *w,
		double *alpha, double *eta, double *tau2, int N, int T)
{
	alpha[0] = alpha0;
	eta[0] = eta0;
	tau2[0] = tau20;
	for (int i = 0; i < N; ++i)
		w[i] = rnorm(2.0, 1.0);
}

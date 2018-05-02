#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "negpotential.h"
#include "double_metropolis.h"
#include "regular_metropolis.h"

/**
 * @file dm_test.c
 * @author Naihui Zhou (nzhou@iastate.edu)
 * to run gcc -Wall -pedantic -o dmtest 
 * 	dm_test.c negpotential.c double_metropolis.c -lm -lRmath
 *
 **/
//typedef double (* pdf)(double, double *);
//pdf functions, fisrt  x, second is pointer to other parameters;
//typedef double (*negp) (double *, double);
//first x, second parameter. only one parameter;
//typedef void (*auxiliary)(int, double *, double *,int **, double *);
//generate auxiliary variable y, given special parameter;


void allocate(double **w, double ** alpha, 
		double ** eta, double ** tau2,int N, int T){
	*w = aligned_alloc(32, T * N * sizeof **w);
	*alpha = aligned_alloc(32, T*sizeof **alpha);
	*eta = aligned_alloc(32, T*sizeof **eta);
	*tau2 = aligned_alloc(32, T*sizeof **tau2);
	if ((*w == NULL)| (*alpha == NULL)| (*eta == NULL) | (*tau2 == NULL)){
		printf("Error allocating memory\n");
		exit(EXIT_FAILURE);
	}
}

void allocate_column(double *w, double **w_bycol, int N, int T){
	for (int i = 0;i<T;i++){
		w_bycol[i] = (double *)w + i*N;
	}
}
void free_memory(double *w, double **w_bycol, double *alpha, double *eta, double *tau2){
	free(w);
	free(alpha);
	free(eta);
	free(tau2);
	free(w_bycol[0]);
	free(w_bycol);
}


/**
 * w is initialized to have rnorm(2,1)
 * w has N rows and T columns
 * w[i,j] = w[i*N+j]
 * we fill out the first column (iteration 0)
 **/
void initialize(double alpha0, double eta0, double tau20, double *w, double *alpha, double *eta, double *tau2, int N, int T){
	alpha[0] = alpha0;
	eta[0] = eta0;
	tau2[0] = tau20;
	for (int i =0;i<N;i++){
		w[i] = rnorm(2,1);
	}
}

/**
 * print to R readable format
 * temporary for this test script
 **/

void print(int T, double* alpha, double* eta, double* tau2){
	printf("alpha=c(");
	for (int i=0;i<T;i++){
		printf("%.4g,",alpha[i]);
	}
	printf(")\n");
	printf("eta=c(");
	for (int i=0;i<T;i++){
		printf("%.4g,",eta[i]);
	}
	printf(")\n");
	printf("tau2=c(");
	for (int i=0;i<T;i++){
		printf("%.4g,",tau2[i]);
	}
	printf(")\n");
}

int main(void){
	int N = 100; //data size
	int T = 2000; //number of iterations;
	//int M = 1; // number of starting value sets;
	//Not doing scale reduction factor;
	//But keeping the option open;
	double * w = NULL;
	double **w_bycol = NULL;
	w_bycol = malloc(T*sizeof **w_bycol);
	double * alpha = NULL;
	double * eta = NULL;
	double * tau2 = NULL;

	
	double *y = NULL; // y should come from R;
	int **neighbor;//neighbor should come from R;
	//use random data for y in this test script;
	//use random data for neighbor in this test script;
	y = malloc(N*sizeof(double));
	neighbor = malloc(N*sizeof(int *));
	for (int i =0; i<N;i++)
		neighbor[i] = malloc(N*sizeof(int));

	if (!y  | !neighbor) {
		printf("Error allocating memory!\n");
		exit(EXIT_FAILURE);
	}
	for (int i = 0; i < N; ++i) {
		y[i] = rnorm(0,0.2);
		//printf("%g", y[i]);
		for (int j = 0; j<N;j++){
			neighbor[i][j] = rbinom(1,0.3);
			//printf("neighbor matrix %d %d is %d", i, j, neighbor[i][j]);
		}
	}
	//Done giving random y and neighbor data;
	
	allocate( &w, &alpha, &eta, &tau2, N, T);
	allocate_column(w, w_bycol, N, T);
	initialize(0.5,-0.16,0.5, w, alpha, eta, tau2, N, T);
	
	//**************************************************;
	//parameters for metropolis
	//vars are the variance for the generated values in metropolis
	//bounds are the upper and lower bounds of the prior uniform
	//distributions for the three parameters
	double vars[4] = {0.1,0.4,0.13,0.7};
	double alpha_bounds[2] = {0.0,10.0};
	double eta_bounds[2] = {-0.199,0.154};
	double tau2_bounds[2] = {0.0,10.0};
	//
	//
	//**************************************************;
	int t; //iteration counter; 
	/*
	within iteration t:
		- step1 Use rm (regular metropolis) to generate posterior w;
		- step2 Use dm (double metropolis) to generate posterior alpha
		- step2 Use dm to generate posterior eta
		- step3 Use dm to generate posterior tau2
	*/
	double new_alpha, new_eta, new_tau2;
	for (t = 0;t<T;t++){
		printf("MC Iteration %d\n", t+1);
		//step1;
		metropolis_for_w_univar(t,N, w_bycol, y, vars[0]);
		//step2 (alpha);
		new_alpha =dm_step1(alpha[t],prior_alpha,vars[1],alpha_bounds);
		alpha[t+1] =dm_step2_t_alpha(t,w_bycol,alpha,eta,tau2,N, T, new_alpha, neighbor);
		//step3 (eta);
		new_eta =dm_step1(eta[t],prior_eta,vars[2],eta_bounds);
		eta[t+1] =dm_step2_t_eta(t,w_bycol,alpha,eta,tau2,N, T, new_eta, neighbor);
		//step4 (tau2);
		new_tau2 =dm_step1(tau2[t],prior_tau2,vars[3],tau2_bounds);
		tau2[t+1] =dm_step2_t_tau2(t,w_bycol,alpha,eta,tau2,N, T, new_tau2, neighbor);
	}
	print(T, alpha, eta, tau2);
	free_memory(w,w_bycol,alpha,eta,tau2);
	free(y);
	free(neighbor);
}



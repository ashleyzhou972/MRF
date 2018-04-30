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
l
typedef double (* pdf)(double, double *);
//pdf functions, fisrt  x, second is pointer to other parameters;
typedef double (*negp) (double *, double);
//first x, second parameter. only one parameter;
typedef void (*auxiliary)(int, double *, double *,int **, double *);
//generate auxiliary variable y, given special parameter;


<<<<<<< HEAD
void allocate(double **w, double **w_bycol, double ** alpha, 
		double ** eta, double ** tau2,int N, int T){
	*w = aligned_alloc(32, T * N * sizeof **w);
	*alpha = aligned_alloc(32, T*sizeof **alpha);
	*eta = aligned_alloc(32, T*sizeof **eta);
	*tau2 = aligned_alloc(32, T*sizeof **tau2);
	if (!(*w & *alpha & *eta & *tau2)){
		printf("Error allocating memory\n");
		exit(EXIT_FAILURE);
	}
	for (int i=0;i<T;i++){
		w_bycol[i] = *w+i*N;
	}
}

void allocate_column(double *w, double **w_bycol, int N, int T){
	for (int i = 0;i<T;i++){
		w_bycol[i] = *w + i*N;
	}
}
void free(double *w, double **w_bycol, double *alpha, double *eta, double *tau2){
=======
void allocate(double ** w, double ** alpha, double ** eta, double ** tau2,
		int N, int T){
	*w = aligned_alloc(64, T * N * sizeof **w);
	*alpha = aligned_alloc(64, T*sizeof **alpha);
	*eta = aligned_alloc(64, T*sizeof **eta);
	*tau2 = aligned_alloc(64, T*sizeof **tau2);
	if (!(w & alpha & eta & tau2)){
		printf("Error allocating memory\n");
		exit(EXIT_FAILURE);
	}
}

void free(double ** w, double **alpha, double ** eta, double ** tau2){
>>>>>>> 031c08d20fc11466614949ba4ab88f6377e47b8b
	free(w);
	free(alpha);
	free(eta);
	free(tau2);
<<<<<<< HEAD
	free(w_bycol[0]);
	free(w_bycol);
=======
>>>>>>> 031c08d20fc11466614949ba4ab88f6377e47b8b
}


/**
 * w is initialized to have rnorm(2,1)
 * w has N rows and T columns
<<<<<<< HEAD
 * w[i,j] = w[i*N+j]
=======
>>>>>>> 031c08d20fc11466614949ba4ab88f6377e47b8b
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

<<<<<<< HEAD

=======
>>>>>>> 031c08d20fc11466614949ba4ab88f6377e47b8b

int main(void){
	int N = 100; //data size
	int T = 2000; //number of iterations;
	int M = 1; // number of starting value sets;
	//Not doing scale reduction factor;
	//But keeping the option open;
	double * w = NULL;
	double **w_bycol = NULL;
	w_bycol = malloc(T*sizeof **w_bycol);
	double * alpha = NULL;
	double * eta = NULL;
	double * tau2 = NULL;
	double y[N]; // y should come from R;
	//use random data for y in this test script;
	allocate( &w, &alpha, &eta, &tau2, N, T);
	allocate_column(w, w_bycol, N, T);
	initialize(0.5,-0.16,0.5, w, alpha, eta, tau2, N, T);
	int t; //iteration counter; 

	/*
	within iteration t:
		- step1 Use rm (regular metropolis) to generate posterior w;
		- step2 Use dm (double metropolis) to generate posterior alpha
		- step2 Use dm to generate posterior eta
		- step3 Use dm to generate posterior tau2
	*/
<<<<<<< HEAD
	int i; //iterates through N
	for (t = 0;t<T;t++){
		printf("MC Iteration %d\n", t+1);
		//step1;
		metropolis_for_w_univar(t, w, y, 0.5);
		//step2 (alpha);
		new_alpha =dm_step1(alpha[t],prior_alpha,0.5,alpha_prior_params);
		alpha[t+1] =dm_step2_alpha_t(t,w_bycol,alpha,eta,tau2,N, T, new_alpha, neighbor);
		//step3 (eta);
		//step4 (tau2);
	}
	free(w,alpha,eta,tau2);
	free(w_bycol);
=======
	for (t = 1;t<T;t++){
		printf("MC Iteration %d\n", t);
		


		
	}
>>>>>>> 031c08d20fc11466614949ba4ab88f6377e47b8b

}



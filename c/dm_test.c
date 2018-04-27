#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <negpotential.h>
#include <double_metropolis.h>

/**
 * @file dm_test.c
 * @author Naihui Zhou (nzhou@iastate.edu)
 * to run gcc -Wall -pedantic -o dmtest 
 * 	dm_test.c negpotential.c double_metropolis.c -lm -lRmath
 *
 **/

typedef double (* pdf)(double, double *);
//pdf functions, fisrt  x, second is pointer to other parameters;
typedef double (*negp) (double *, double);
//first x, second parameter. only one parameter;
typedef void (*auxiliary)(int, double *, double *,int **, double *);
//generate auxiliary variable y, given special parameter;


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
	free(w);
	free(alpha);
	free(eta);
	free(tau2);
}


/**
 * w is initialized to have rnorm(2,1)
 * w has N rows and T columns
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


int main(void){
	int N = 100; //data size
	int T = 2000; //number of iterations;
	int M = 1; // number of starting value sets;
	//Not doing scale reduction factor;
	//But keeping the option open;
	double ** w = NULL;
	double * alpha = NULL;
	double * eta = NULL;
	double * tau2 = NULL;
	int t; //iteration counter; 
	/*
	within iteration t:
		- Use rm (regular metropolis) to generate posterior w;
		- Use dm (double metropolis) to generate posterior alpha
		- Use dm to generate posterior eta
		- Use dm to generate posterior tau2
		- update w, alpha, eta, tau2
	*/
	for (t = 1;t<T;t++){
		printf("MC Iteration %d\n", t);
		


		
	}

}


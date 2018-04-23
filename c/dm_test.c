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

}


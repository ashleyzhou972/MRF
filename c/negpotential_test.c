#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "negpotential.h"
#include <time.h>
/**
 * @file negpotential_test.c
 * @author Naihui Zhou (nzhou@iastate.edu)
 * to run gcc -Wall -pedantic -o negtest 
 * 	negpotential_test.c negpotential.c -lm -lRmath
 **/



int main(void)
{
	int N = 20000;
	double *y = NULL;
	double *ystar = NULL;
	srand(time(NULL));
	unsigned int seed1 = rand();
	unsigned int seed2 = rand();
	int **neighbor = NULL;

	set_seed(seed1, seed2);

	int i;

	y = malloc(N*sizeof(double));
	ystar = malloc(N*sizeof(double));
	neighbor = malloc(N*sizeof(int *));

	for (i =0; i<N;i++)
		neighbor[i] = malloc(N*sizeof(int));

	if (!y & !ystar & !neighbor) {
		printf("Error allocating memory!\n");
		exit(EXIT_FAILURE);
	}

	for (i = 0; i < N; ++i) {
		y[i] = rnorm(0,1.5);
		ystar[i] = 0.0;
		//printf("%g", y[i]);
		for (int j = 0; j<N;j++){
			neighbor[i][j] = rbinom(1,0.3);
			//printf("neighbor matrix %d %d is %d", i, j, neighbor[i][j]);
		}
	}

	double out;

	out = negpotential(y, ystar, N, neighbor, 2.0, 0.15,5.0);
	printf("result is %g\n", out);
	free(y);
	free(ystar);
	for (i = 0; i < N; ++i)
		free(neighbor[i]);
	free(neighbor);
	return 0;
}

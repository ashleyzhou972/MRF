#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "../negpotential.h"
#include <time.h>
/**
 * @file negpotential_test.c
 * @author Naihui Zhou (nzhou@iastate.edu)
 * to run gcc -Wall -pedantic -o negtest negpotential_test.c ../negpotential.c -lm -lRmath -fopenmp -DOPENMP
 **/

int main(void)
{
	int N = 100;
	double *y = NULL;
	double *ystar = NULL;
	unsigned int seed1 = rand();
	unsigned int seed2 = rand();
	int **neighbor = NULL;
	int i;

	srand(time(NULL));
	set_seed(seed1, seed2);

	y = malloc(N*sizeof(double));
	ystar = malloc(N*sizeof(double));
	neighbor = malloc(N*sizeof(int *));

	for (i = 0; i < N; ++i)
		neighbor[i] = malloc(N*sizeof(int));

	if (!y | !ystar | !neighbor) {
		printf("Error allocating memory!\n");
		exit(EXIT_FAILURE);
	}

	for (i = 0; i < N; ++i) {
		y[i] = rnorm(0.0, 0.2);
		ystar[i] = 0.0;
		for (int j = 0; j < N; ++j)
			neighbor[i][j] = rbinom(1, 0.3);
	}

	double out;

	out = negpotential(y, ystar, N, neighbor, 2.0, 0.15, 5.0);
	printf("result is %g\n", out);
	free(y);
	free(ystar);
	for (i = 0; i < N; ++i)
		free(neighbor[i]);
	free(neighbor);
	return 0;
}

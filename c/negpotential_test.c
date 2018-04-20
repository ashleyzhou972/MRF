#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <negpotential.h>
/**
 * @file negpotential_test.c
 * @author Naihui Zhou (nzhou@iastate.edu)
 * to run gcc -Wall -pedantic -o negtest 
 * 	negpotential_test.c negpotential.c -lm -lRmath
 *
 **/



int main(void){
	double * y = NULL;
	double * ystar = NULL;
	unsigned int seed1 = rand();
	unsigned int seed2 = rand();
	int ** neighbor = NULL;
	set_seed(seed1, seed2);
	int i;
	y = malloc(100*sizeof(double));
	ystar = malloc(100*sizeof(double));
	neighbor = malloc(100*sizeof(int *));
	for (i =0; i<100;i++) neighbor[i] = malloc(100*sizeof(int));
	if (!y & !ystar & !neighbor)  {
		
		printf("Error allocating memory!\n");
		exit(EXIT_FAILURE);
	} 

	for (i = 0; i<100;i++){
		y[i] = rnorm(0,1.5);
		ystar[i] = 0.0;
		//printf("%g", y[i]);
		for (int j = 0; j<100;j++){
			neighbor[i][j] = rbinom(1,0.5);
			//printf("neighbor matrix %d %d is %d", i, j, neighbor[i][j]);
		}
	}
	double out;
	out = negpotential(y, ystar, 100, neighbor, 2.0, 0.15,5.0);
	printf("result is %g\n", out);
	free(y);
	free(ystar);
	for (i=0;i<100;i++) free(neighbor[i]);
	free(neighbor);
	return 0;
}

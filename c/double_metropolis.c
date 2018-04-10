/**
 * 
 * Initial attempt at the double metropolis algorithm 
 * @author Naihui Zhou {nzhou@iastate.edu}
 * to run gcc -Wall -pedantic -o dm double_metropolis.c -lm -lRmath
 **/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

/**
 * This is the Q function
 **/
double negpotential(double* w, double alpha, double eta, double tau2,); 
double Hij(int i, int j, double alpha, double eta, double tau2);


int main(void){
	return 0;
}

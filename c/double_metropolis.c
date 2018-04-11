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

double negpotential(double* w, double alpha, double eta, double tau2,); 
double Hij(double y_i, double  ystar_i, double y_j, double ystar_j, 
		int i, int j, double alpha, double eta, double tau2);

double data_pdf();
double mixing_pdf();
double prior_pdf();

double posterior_pdf_w();
/**
 * Each of the parameters alpha, eta and tau2
 * has two types of posteriors:
 * 	- The multivariate normal -1
 * 	- The negpotential function -2
 * For the negpotential function type posterior,
 * there are two methods to compute the intractable constant
 * 	- The double- metropolis algorithm
 * 	- importance sampling -3
 * Only case 2 will be in this script
 **/
double posterior_pdf_alpha();
double posterior_pdf_eta();
double posterior_pdf_tau2();

void dm_step1();
void dm_step2();
void dm_step3();



int main(void){
	return 0;
}

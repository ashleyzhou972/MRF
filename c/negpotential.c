#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

/**
 * @file negpotential.c
 * @author Naihui Zhou (nzhou@iastate.edu)
 * Function to calculate negpotential function Q
 *
 **/

double negpotential(double * y, double * ystar, int size_y, double ** neighbor,
	       	double alpha, double eta, double tau2);
double H_i(double * y, int i,  double * ystar, int size_y, double ** neighbor, 
		double alpha, double eta, double tau2);
double H_ij(double * y, int i, int j,  double * ystar, int size_y ,
		double ** neighbor, double alpha, double eta, double tau2);
double vector_multiplication(double * array1, double * array2, int size);


double negpotential(double * y, double * ystar, int size_y, double ** neighbor,
	       	double alpha, double eta, double tau2){
	double summand_i = 0;
	double summand_ij = 0;
	double result;
	for (int i = 0; i < size_y; i++){
		summand_i += H_i(y, i, ystar, size_y, neighbor, alpha, eta,tau2);
		for (int j = 0; j < size_y; j++){
			summand_ij += H_ij(y, i, j, ystar, size_y, neighbor, 
					alpha, eta, tau2);
		}
	}
	result = summand_i + summand_ij;
	return result;
}




double H_i(double * y, int i,  double * ystar, int size_y, double ** neighbor, 
		double alpha, double eta, double tau2){
	double mu, log_dnorm;
	mu = alpha + eta*vector_multiplication(neighbor[i],ystar-alpha,size_y);
	log_dnorm = dnorm(y[i], mu, tau2, 1)-dnorm(ystar[i],mu,tau2,1);
       		//1 for log = TRUE
	return log_dnorm;

}


double H_ij(double * y, int i, int j,  double * ystar, int size_y ,
		double ** neighbor, double alpha, double eta, double tau2){
	double mu1, mu2 log_dnorm;
	double * ystar_temp = ystar;
	ystar_temp[j] = y[j];
	mu1 = alpha + eta*vector_multiplication(neighbor[i], ystar_temp, size_y);
	mu2 = alpha + eta*vector_multiplication(neighbor[i], ystar,size_y);
	log_dnorm = dnorm(y[i],mu1,tau2,1) + dnorm(ystar[i],mu2,tau2,1)
		- dnorm(ystar[i],mu1,tau2,1) - dnorm(y[i],mu2,tau2,1);
	return log_dnorm;
}


double vector_multiplication(double * array1, double * array2, int size){
	double summand = 0;
	for (int i=0;i<size;i++){
		summand += array1[i]*array2[i];
	}
	return summand;
}

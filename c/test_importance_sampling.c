#include <stdio.h>
#include <stdlib.h>

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "importance_sampling.h"

typedef double(*func)(double);


// Test P(Y>=0.3) given Y~N(0,1) truth == 0.00139

double h(double Y)
{
	fprintf(stdout,"In h\n");
	if (Y > 3.0)
		return 1.0;
	else
		return 0.0;
}

double g(double Y)
{
	fprintf(stdout,"In g\n");
	return dnorm(Y,4.0,1.0,0);
}

double G(double X)
{
	fprintf(stdout,"In G\n");
	return rnorm(4.0,1.0);
}

double f(double Y)
{
	fprintf(stdout,"In f\n");
	return dnorm(Y,0.0,1.0,0);
}


int main(void)
{
	double E = 0.0;
	int const N = 1000;
	
	double temp = h(0.0);
	temp = g(0.0);
	temp = G(0.0);
	temp = f(0.0);
	
	func *ptr2h = &h;
	func *ptr2g = &g;
	func *ptr2G = &G;
	func *ptr2f = &f;

	E = do_important_sampling(ptr2h, ptr2f, ptr2g, ptr2G, N); // NOT WORKING!!!
/*	for(int k = 0; k < N; ++k) {*/
/*		double X = runif(0.0,1.0);*/
/*		double Y = G(X);*/
/*		E += h(Y)*f(Y)/g(Y);*/
/*	}*/
/*	E /= N;*/
	fprintf(stdout,"E = %g\n",E);
	return 0;
}

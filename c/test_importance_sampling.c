#include <stdio.h>
#include <stdlib.h>

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "importance_sampling.h"

typedef double(*func)(double);


// Test P(Y>=0.3) given Y~N(0,1) truth == 0.00139

double h(double Y)
{
	if (Y > 3.0)
		return 1.0;
	else
		return 0.0;
}

double g(double Y)
{
	return dnorm(Y,4.0,1.0,0);
}

double G(double X)
{
	return rnorm(4.0,1.0);
}

double f(double Y)
{
	return dnorm(Y,0.0,1.0,0);
}


int main(void)
{
	int flag;
	double E = 0.0;
	double Var = 0.0;
	double Output[2];
	int const N = 100000;

	func ptr2h = &h;
	func ptr2g = &g;
	func ptr2G = &G;
	func ptr2f = &f;

	flag = do_important_sampling(ptr2h, ptr2f, ptr2g, ptr2G, N, Output);
	if (flag == 0) {
		E = Output[0];
		Var = Output[1];
		fprintf(stdout,"E = %g\n",E);
		fprintf(stdout,"Var = %g\n",Var);
	} else
		fprintf(stdout,"Something Bad Happened...\n");
	return 0;
}

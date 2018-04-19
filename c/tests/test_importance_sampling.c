#include <stdio.h>
#include <stdlib.h>

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "../importance_sampling.h"

// typedef double(*func)(double);

double h(double Y)
{
	if (Y > 3.0)
		return 1.0;
	else
		return 0.0;
}

// Known probability function that we can sample from
double g(double Y) {return dnorm(Y, 4.0, 1.0, 0); }

// CDF of g that we can sample from quickly...
double G(double X)
{
	// X ~ U(0,1)
	++X; // ++X to remove warning...
	return rnorm(4.0, 1.0);
}

// Density we care about
double f(double Y) {return dnorm(Y, 0.0, 1.0, 0); }


int main(void)
{
	int flag;
	double Output[2];
	int const N = 1000;

	fprintf(stdout,"Test P(Y>=3.0) given Y~N(0,1)\n");

	flag = do_important_sampling(h, f, g, G, N, Output);
	if (flag == 0) {
		double E = Output[0];
		double Var = Output[1];
		double Tru_V = pnorm(3.0, 0.0, 1.0, 0, 0);
		double Err_V = fabs(E - Tru_V);

		fprintf(stdout, "E = %g\n", E);
		fprintf(stdout, "Var = %g\n", Var);
		fprintf(stdout, "Truth = %g\n", Tru_V);
		fprintf(stdout, "Error = %g\n", Err_V);
	} else
		fprintf(stdout, "Something Bad Happened...\n");
	return 0;
}

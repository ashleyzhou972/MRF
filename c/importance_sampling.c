// Need to use random tool box from R
#include <stdlib.h>
#include <stdio.h> // For debuging...

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "importance_sampling.h"

double do_important_sampling(func *h, func *f, func *g, func *G, int numsamples)
{
	double E = 0.0;

	for(int k = 0; k < numsamples; ++k) {
		fprintf(stdout,"E = %g\n",E);
		double X = runif(0.0,1.0);
		fprintf(stdout,"X = %g\n",X);
		double Y = (*G)(X);
		fprintf(stdout,"Y = %g\n",Y);
		E += (*h)(Y) * (*f)(Y) / (*g)(Y);
	}

	E /= numsamples;

	return E;
}


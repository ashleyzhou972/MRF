// Need to use random tool box from R
#include <stdlib.h>
#include <stdio.h> // For debuging...

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "importance_sampling.h"

int do_important_sampling(func h, func f, func g, func G, int numsamples, double *output)
{
	double X, Y, Enp1;
	double E = 0.0;
	double E_old = 0.0;
	double V = 0.0;
	double V_old = 0.0;

	for(int k = 0; k < numsamples; ++k) {
		E_old = E;
		V_old = V;
		X = runif(0.0,1.0);
		Y = (*G)(X);

		Enp1 = (*h)(Y) * (*f)(Y) / (*g)(Y);
		E = k*E_old + Enp1;
		E /= k+1;

		V = (E_old - Enp1) * (E_old - Enp1);
		V /= k+1;
		V += V_old;
		V *= (double)k / (k + 1);
	}

/*	E /= numsamples;*/
	output[0] = E;
	output[1] = V;

	return 0;
}


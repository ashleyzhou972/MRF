#include <stdio.h>
#include <stdlib.h>

#include <Rmath.h>

int main(void) {
	double rand;
	rand = rnorm(0,1);
	printf("N(0,1) random number is %g\n", rand);
	return 0;
}

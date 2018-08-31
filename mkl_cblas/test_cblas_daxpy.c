#include <stdio.h>
#include <stdlib.h>

#include <mkl.h>

void scalar_vector_multiplication(int n, double scalar, double *out, double *vector)
{
	int incx = 1.0;
	int incy = 1.0;
	cblas_daxpy(n, scalar, vector, incx, out, incy);
}

/**
 * scalar repeated n times, to form a identical vector of a;
 * plus another input vector
 * @param vector is updated;
 **/
void scalar_vector_summation(int n, double scalar, double *vector)
{

	double x[n];
	for (int i = 0; i<n ; ++i) 
		x[i] = 1.0;
	int incx = 1;
	int incy = 1;
	cblas_daxpy(n, scalar, x, incx, vector, incy);
} 
int main(void) {
	int n = 5;
	double *out;
	double *vector;
	out = malloc(n*sizeof(double));
	vector = malloc(n*sizeof(double));
	if (out == NULL || vector == NULL) {
		printf("error allocating memory\n");
	}
	else {
		vector[0] =1.0;
		vector[1] =1.0;
		vector[2] =2.0;
		vector[3] =3.0;
		vector[4] =4.0;
	}
	//print matrix and vecotr;
	printf("The input vector is \n");
	for (int i = 0; i<n; ++i) {
		printf("%.5f ", vector[i]);
	}	
	printf("\n");
	double scalar = 2.0;
	scalar_vector_multiplication(n, scalar, out, vector);
	scalar_vector_summation(n, 5, out);
	free(vector);
	printf("The output vector is \n");
	for (int i =0; i<n;++i) {
		printf("%.5g\n", out[i]);
	}
	printf("\n");
	free(out);
}

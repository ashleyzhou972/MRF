#include <stdio.h>
#include <stdlib.h>

#include <mkl.h>
#include <Rmath.h>

double vector_dot_product(int n, int *vector1, double *vector2) {
	int incx = 1;
	int incy = 1;
	//convert int vector to double;
	double *dvector1;
	dvector1 = malloc(n*sizeof(double));
	for (int i = 0;i<n;++i) {
		dvector1[i] = (double) vector1[i];
	}
	double res;
	res = cblas_ddot(n, dvector1, incx, vector2, incy);
	free(dvector1);
	return res;
}
int main(void) {
	int n = 5;
	int *vector1;
	double *vector2;
	vector1 = malloc(n*sizeof(int));
	vector2 = malloc(n*sizeof(double));
	if (vector1 == NULL || vector2 == NULL) {
		printf("error allocating memory\n");
	}
	else {
		vector1[0] =1;
		vector1[1] =1;
		vector1[2] =2;
		vector1[3] =3;
		vector1[4] =4;
		vector2[0] =1.0;
		vector2[1] =1.0;
		vector2[2] =2.0;
		vector2[3] =3.0;
		vector2[4] =4.0;

	}
	//print matrix and vecotr;
	printf("The input vector 1 is \n");
	for (int i = 0; i<n; ++i) {
		printf("%d ", vector1[i]);
	}	
	printf("\n");
	printf("The input vector 2 is \n");
	for (int i = 0; i<n; ++i) {
		printf("%.5f ", vector2[i]);
	}printf("\n");
	double out;
	out = vector_dot_product(n, vector1, vector2);
	printf("The output value is \n");
	printf("%.5g\n", out);
	printf("\n");
	free(vector1);
	free(vector2);
	double rand;
	rand = rnorm(0,1);
	printf("The random N(0,1) number is %.5f", rand);
}

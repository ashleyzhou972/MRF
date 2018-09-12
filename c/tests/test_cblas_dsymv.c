#include <stdio.h>
#include <stdlib.h>

#include <mkl.h>

void matrix_vector_multiplication(int n, double *out, int *matrix, double *vector)
{
      	MKL_INT         lda, ldb, ldc;
      	MKL_INT         rmaxa, cmaxa, rmaxb, cmaxb, rmaxc, cmaxc;
      	double          alpha, beta;
      	double         *a, *b, *c;
      	CBLAS_LAYOUT    layout;
      	CBLAS_SIDE      side;
      	CBLAS_UPLO      uplo;
      	MKL_INT         ma, na, typeA;
	alpha = 1.0;
	beta = 0.0;
	layout = CblasColMajor;	
	side = CblasLeft;
	uplo = CblasLower;
	layout = CblasColMajor;
	int incx = 1.0;
	int incy = 1.0;
	double *dmatrix = (double *) matrix;
	cblas_dsymv(layout, uplo, n, alpha, dmatrix, lda, vector, incx, beta, out, incy);
}

int main(void) {
	int n = 5;
	double *out;
	double *vector;
	int *matrix;
	out = malloc(n*sizeof(double));
	vector = malloc(n*sizeof(double));
	matrix = malloc(n*n*sizeof(int));
	if (out == NULL || vector == NULL || matrix == NULL) {
		printf("error allocating memory\n");
	}
	else {
		matrix[0] = 0.0;
		matrix[5] = 1.0;
		matrix[6] = 0.0;
		matrix[10] = 0.0;
		matrix[11] = 1.0;
		matrix[12] = 0.0;
		matrix[15] = 1.0;
		matrix[16] = 0.0;
		matrix[17] = 0.0;
		matrix[18] = 0.0;
		matrix[20] = 1.0;
		matrix[21] = 1.0;
		matrix[22] = 1.0;
		matrix[23] = 0.0;
		matrix[24] = 0.0;
		vector[0] =1.0;
		vector[1] =1.0;
		vector[2] =2.0;
		vector[3] =3.0;
		vector[4] =4.0;
	}
	matrix_vector_multiplication(n, out, matrix, vector);
	for (int i; i<n;++i) {
		printf("%g\n", out[i]);
	}
	free(out);
	free(vector);
	free(matrix);
}

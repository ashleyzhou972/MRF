#include <stdio.h>
#include <stdlib.h>

#include <mkl.h>

void matrix_vector_multiplication(int n, double *out, int *matrix, double *vector)
{
      	MKL_INT         lda;
      	double          alpha, beta;
      	CBLAS_LAYOUT    layout;
      	CBLAS_SIDE      side;
      	CBLAS_UPLO      uplo;
	alpha = 1.0;
	beta = 0.0;
	layout = CblasRowMajor;	
	uplo = CblasLower;
	lda = n ;
	int incx = 1.0;
	int incy = 1.0;
	//cast input matrix to double ;
	double * dmatrix;
	dmatrix = malloc(n*n*sizeof(double));
	for (int i = 0; i<n*n; ++i) {
		dmatrix[i] = (double) matrix[i];
	}
	cblas_dsymv(layout, uplo, n, alpha, dmatrix, lda, vector, incx, beta, out, incy);
	free(dmatrix);
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
		matrix[0] = 0;
		matrix[1] = 1;
		matrix[2] = 0;
		matrix[3] = 1;
		matrix[4] = 1;
		matrix[5] = 1;
		matrix[6] = 0;
		matrix[7] = 1;
		matrix[8] = 0;
		matrix[9] = 1;
		matrix[10] = 0;
		matrix[11] = 1;
		matrix[12] = 0;
		matrix[13] = 0;
		matrix[14] = 1;
		matrix[15] = 1;
		matrix[16] = 0;
		matrix[17] = 0;
		matrix[18] = 0;
		matrix[19] = 0;
		matrix[20] = 1;
		matrix[21] = 1;
		matrix[22] = 1;
		matrix[23] = 0;
		matrix[24] = 0;
		vector[0] =1.0;
		vector[1] =1.0;
		vector[2] =2.0;
		vector[3] =3.0;
		vector[4] =4.0;
	}
	//print matrix and vecotr;
	printf("The input int matrix is \n");
	for (int i=0;i<n;++i) {
		for (int j=0;j<n;++j) {
			printf("%d ", matrix[i*n + j]);
		}
		printf("\n");
	}
	//convert int matrix to double via loop;
	double *dmatrix;
	dmatrix = malloc(n*n*sizeof(double));
	for (int i=0;i<n*n;++i) {
		dmatrix[i] = (double) matrix[i];
	}

	printf("The input double matrix is \n");
	for (int i=0;i<n*n;++i) {
		printf("%.5g ", dmatrix[i]);
	}
	printf("\n");

	printf("The input vector is \n");
	for (int i = 0; i<n; ++i) {
		printf("%.5f ", vector[i]);
	}	
	printf("\n");
	matrix_vector_multiplication(n, out, matrix, vector);
	free(vector);
	free(matrix);
	free(dmatrix);
	printf("The output vector is \n");
	for (int i =0; i<n;++i) {
		printf("%.5g\n", out[i]);
	}
	printf("\n");
	//printf("hello world\n");
	
	free(out);
	}

#include <stdio.h>
#include "mkl.h"

void cblas_dsymm (const CBLAS_LAYOUT Layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const MKL_INT m, const MKL_INT n, const double alpha, const double *a, const MKL_INT lda, const double *b, const MKL_INT ldb, const double beta, double *c, const MKL_INT ldc);


void matrix_vector_multiplication(double **matrix, double **vector, double **output, int N) {
	double alpha = 1.0;
	double beta = 0.0;
	
	cblass_dsymm(
}

int main(void){
	double *A =  

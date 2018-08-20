#include <stdio.h>
#include <stdlib.h>

#include "mkl.h"

#include "mkl_cblas.h"
#include "mkl_example.h"

int main(void) {
	MKL_INT         m, n;
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
	layout = CblasRowMajor;	
	side = CblasLeft;
	uplo = CblasLower;
	m = 5;
	n = 1;
	if( side == CblasLeft ) {
		rmaxa = m + 1;
	  	cmaxa = m;
	  	ma    = m;
	  	na    = m;
      	} else {
	  	rmaxa = n + 1;
	  	cmaxa = n;
	  	ma    = n;
	  	na    = n;
      }
      	rmaxb = m + 1;
      	cmaxb = n;
      	rmaxc = m + 1;
      	cmaxc = n;
      	a = (double *)mkl_calloc(rmaxa*cmaxa, sizeof( double ), 64);
      	b = (double *)mkl_calloc(rmaxb*cmaxb, sizeof( double ), 64);
      	c = (double *)mkl_calloc(rmaxc*cmaxc, sizeof( double ), 64);
	/**
 * 	input a values here
  	**/
	//a = {1.0, 0.0, 1.0,0.0, 0.0, 1.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0};
	a[0] = 0.0;
	a[5] = 1.0;
	a[6] = 0.0;
	a[10] = 0.0;
	a[11] = 1.0;
	a[12] = 0.0;
	a[15] = 1.0;
	a[16] = 0.0;
	a[17] = 0.0;
	a[18] = 0.0;
	a[20] = 1.0;
	a[21] = 1.0;
	a[22] = 1.0;
	a[23] = 0.0;
	a[24] = 0.0;
	b[0] =1.0;
	b[1] =1.0;
	b[2] =2.0;
	b[3] =3.0;
	b[4] =4.0;
	//b = {1.0,1.0,2.0,3.0,4.0};
      	if( a == NULL || b == NULL || c == NULL ) {
	  	printf( "\n Can't allocate memory for arrays\n");
	  	return 1;
      	}
      	if( layout == CblasRowMajor ) {
	 	lda=cmaxa;
	 	ldb=cmaxb;
	 	ldc=cmaxc;
      	} else {
	 	lda=rmaxa;
	 	ldb=rmaxb;
	 	ldc=rmaxc;
      	}
      	cblas_dsymm(layout, side, uplo, m, n, alpha, a, lda,
                  b, ldb, beta, c, ldc);
/**
 * print matrices;
 **/
	typeA = LOWER_MATRIX;
      	PrintArrayD(&layout, FULLPRINT, typeA, &ma, &na, a, &lda, "A");
     	PrintArrayD(&layout, FULLPRINT, GENERAL_MATRIX, &m, &n, b, &ldb, "B");
      	PrintArrayD(&layout, FULLPRINT, GENERAL_MATRIX, &m, &n, c, &ldc, "C");
	printf("\n");
	mkl_free(a);
	mkl_free(b);
	mkl_free(c);
}




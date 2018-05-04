#include <Rinternals.h>


SEXP double_metropolis(SEXP T_in, SEXP y_in, SEXP neighbor_in, SEXP vars_in, SEXP bounds_alpha, SEXP bounds_eta, SEXP bounds_tau2){
	//declaration of variables;
	int N,T;
	int *T_p;
	int *y_int;
	//double *y;
	int *neighbor_1d;
	int dim_neighbor;
	double *vars;
	double *b_alpha, *b_eta, *b_tau2;
	//Check incoming SEXP types and extract data
	if (!isInteger(T_in))
		error("[ERROR] First argument must be interger");
	
	else {
		T_p = INTEGER(T_in);
		T = *T_p;
	}
	if (!isInteger(y_in) || !isVector(y_in))
		error("[ERROR] Second argument must be integer vector");
	else {
		y_int = INTEGER(y_in);
		N = length(y_in);
	}
	int neighbor_2d[N][N];
	if (!isInteger(neighbor_in) || !isVector(neighbor_in))
		error("[ERROR] Third argument must be integer vector");
	else {
		neighbor_1d = INTEGER(neighbor_in);
		dim_neighbor = length(neighbor_in);
		if (dim_neighbor!=N*N)
			error("[ERROR] Size of neighborhood matrix does not agree with size of data vector");
		else {
			for (int i = 0;i<N;i++) {
				for (int j = 0;j<N;j++) 
					neighbor_2d[i][j] = neighbor_1d[j*N+i];
			}
		}
	}
	if (!isReal(vars_in) || !isVector(vars_in) || !length(vars_in)==4)
		error("[ERROR] Fourth argument must be real vector with 4 elements");
	else 
		vars = REAL(vars_in);
	if (!isReal(bounds_alpha) || !isReal(bounds_eta) || !isReal(bounds_tau2) 
			|| !isVector(bounds_alpha) || !isVector(bounds_eta) ||
			!isVector(bounds_tau2)) 
		error("[ERROR] last three arguments must be real vectors\n");
	else {
		if (!length(bounds_alpha)==2 || !length(bounds_eta)==2 ||
					!length(bounds_tau2)==2)
			error("[ERROR] last three arguments must each be of length 2\n");

		else {
			b_alpha = REAL(bounds_alpha);
			b_eta = REAL(bounds_eta);
			b_tau2 = REAL(bounds_tau2);
		}
	}
	Rprintf("size of data vector is %d\n" , N);
	Rprintf("number of iterations is %d\n" , T);
	/*
	Rprintf("data vector is :\n");
	for (int i =0;i<N;i++) {
		Rprintf("%d " , y_int[i]);
	}
	Rprintf("\n");
	Rprintf("first column of neighbor matrix is :\n");
	for (int i=0;i<N;i++){
		Rprintf("%d ",neighbor_2d[i][0]);
	}
	Rprintf("\n");
	*/
	return R_NilValue;
}


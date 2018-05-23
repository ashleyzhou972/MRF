#ifndef REGULAR_METROPOLIS_H
#define REGULAR_METROPOLIS_H

/**
 * @file regular_metropols.h
 * @author Naihui Zhou (nzhou@iastate.edu)
 *
 **/
typedef double (*pdf)(double, double *);
//pdf functions, fisrt  x, second is pointer to other parameters;



/**
 *
 * For regular metropolis-hastings within Gibbs algorithm
 **/
double r_random_walk_chain(double current_y, double var);
double d_random_walk_chain(double current_y, double proposed_y, double var);
double jump_probability(double current, double proposed, pdf targeti, double *param);

/**
 * Array helper functions
 **/
double vector_multiplication(int *array1, double *array2, int size);

/**
 * regular metropolis for w
 **/

double log_data_density_univar(double y, double w);
double log_mrf_density_univar(int size_w, int i, double w_in, double *w, int **neighbor, double alpha, double eta, double tau2);
double log_sum_density_univar(int size_w, int i, double y, double w_in, double *w, int **neighbor, double alpha, double eta, double tau2);
int metropolis_for_w_univar(int t, int N, double **w, double *y, double var, int **neighbor, double alpha, double eta, double tau2);

#endif

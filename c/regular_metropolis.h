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
 * regular metropolis for w
 **/

double log_data_density_univar(double y, double w);
void metropolis_for_w_univar(int t, int N, double **w, double *y, double var);
double pdf_lddu(double w, double *y);

#endif

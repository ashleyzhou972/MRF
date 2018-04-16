#ifndef IMPORTANCE_SAMPLING
#define IMPORTANCE_SAMPLING

// Expects h, f, g, and G, to be of typedef double (*func) (double);
typedef double(*func)(double);
double do_important_sampling(func *h, func *f, func *g, func *G, int numsamples);

#endif

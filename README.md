# Double Metropolis Algorithm in C 
Naihui Zhou and Matt Schramm
## Goal
The overall goal of the project is using Bayesian inference for a Markov Random Field model fitted to gene expression data on the human genome, where the spatial neighborhood information of each gene is inferred from HiC heatmap.

## Model Specification
See attached model_specification_3.pdf.

## Double Metropolis Algorithm
See attached paper by Liang et al.

## Usage
The functions can be accessed through the R script `dm_call.R` in the `c` folder.
Parameter such as of number MCMC iterations, initial values, ... etc can be changed in this R script. 

For direct manipulation through C, use `dm_test.c` in the `c/test/` folder. (Note that this means using randomly generated data)

Observed data are stored in two `.RData` files.

Note: \
In order to take advantage of OpenMP parallel processing, instead of using 
```
R CMD SHLIB dm_call.c double_metropolis.c negpotential.c regular_metropolis.c
```
to compile the shared R library, use
``` 
 gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG      -fpic  -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c dm_call.c -o dm_call.o -fopenmp -DOPENMP
 gcc -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions -Wl,-z,relro -o dm_call.so dm_call.o double_metropolis.o regular_metropolis.o negpotential.o -lm -lRmath -L/usr/lib/R/lib -lR -lRmath
```
 
 


# Double Metropolis Algorithm in C 
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

Note:    
OpenMP parallel processing is built in. To take advantage of this feature, when compiling the shared R library using 
```
R CMD SHLIB dm_call.c double_metropolis.c negpotential.c regular_metropolis.c
```
instead, use the following,
``` 
gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG      -fpic  -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c dm_call.c -o dm_call.o
gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG      -fpic  -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c regular_metropolis.c -o regular_metropolis.o
gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG      -fpic  -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c double_metropolis.c -o double_metropolis.o
gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG      -fpic  -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c negpotential.c -o negpotential.o -fopenmp -DOPENMP
gcc -std=gnu99 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions -Wl,-z,relro -o dm_call.so dm_call.o regular_metropolis.o double_metropolis.o negpotential.o -L/usr/lib/R/lib -lR -lRmath -fopenmp -DOPENMP
```
 
## Acknowledgment
Matt Schramm contributed to C coding, as part of a course project.
For original course project, see https://git.linux.iastate.edu/nzhou/stat580FinalProject 


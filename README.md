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
Through the `.Call()` function, data are passed to C to process.

For direct manipulation through C, use `dm_test.c` in the `c/test/` folder. (Note that this means using randomly generated data)

Observed data are stored in two `.RData` files.

Note: \\
 


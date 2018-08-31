###Different simulated values and different estimation results


setwd('./')
source('./cblas_dm_call.R')
#source('../../functions_for_metropolis.R')
#source('../../scaleReductionFactor.R')
#source('./simulation_1000.R')
#source('./pseudolikelihood.R')
#Total number of iterations
total_iter = 100
B = 2000
#This observed data are simulated gaussian MRF instead of Poisson
#load observed data
#eta can range from -0.133 to 0.116
###true values:


truevalues = matrix(c(2,-0.1,4,3,0.1,2,1,0.115,7,1,-0.13,2,2,0,5), nrow = 5, byrow = T)
load('../R/simulated_neighbors_1000_poisson.RData')
#conter = 1
for (counter in 1:1){
  alpha_true= truevalues[counter,1]
  eta_true = truevalues[counter,2]
  tau2_true = truevalues[counter,3]
  filename = paste("../R/simulated_1000poisson_",alpha_true,"_", eta_true, "_", tau2_true,".RData",sep = "")
  load(filename)
  cat("Load simulated y complete\n")
  
  N = length(y)
  nb_mat_int = as.integer(sub_neighbor)
  #There are four steps of metropolis in each iteration (including double metropolis)
  #When simulating a new value in each of these metropolis steps, variance of the random walk chain can be specified
  vars = c(0.8,0.3,0.5,0.3)
  #first var is ignored in the gaussian case
  #parameters for the prior distributions (uniform), for alpha, eta and tau^2
  bounds_a = c(0,10)
  bounds_e = c(-0.13,0.11)
  bounds_t = c(0,10)
  #inital guess for alpha, beta and tau^2
  
  inis1 = c(5,0.0,0.5)
  set.seed(34)
  inis2 = rnorm(N, 2, 1)
  ret1 <- dm_call_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis1, inis2)
  jump_count = get_jump_frequency(ret1, total_iter, N)
  print(jump_count)
}



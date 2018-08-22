###Different simulated values and different estimation results


setwd('./')
source('../c/dm_call_cont.R')
source('./functions_for_metropolis.R')
#source('../../scaleReductionFactor.R')
source('./simulation_1000.R')
#source('./pseudolikelihood.R')
#Total number of iterations
total_iter = 100
#B = 2000
#This observed data are simulated gaussian MRF instead of Poisson
#load observed data
#eta can range from -0.133 to 0.116
###true values:


truevalues = matrix(c(2,-0.1,4,3,0.1,2,1,0.115,7,1,-0.13,2,2,0,5), nrow = 5, byrow = T)
load('./simulated_neighbors_1000_poisson.RData')
net = graph_from_adjacency_matrix(sub_neighbor, mode = "undirected")
LOAD_y = TRUE
LOAD_inis = TRUE
date = "20180822"
#conter = 1
for (counter in 1:1){
  alpha_true= truevalues[counter,1]
  eta_true = truevalues[counter,2]
  tau2_true = truevalues[counter,3]
  filename = paste("./simulated_1000poisson_",alpha_true,"_", eta_true, "_", tau2_true,".RData",sep = "")
  if (LOAD_y){
    load(filename)
    cat("Load simulated y complete\n")
  }else{
    y=simulate_y_poisson(net,alpha_true,eta_true,tau2_true,10000)
    save(y,file=filename)
    cat("Simulation done and saved\n")
  }
  
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
  if (!LOAD_inis) {
    inis1 = c(5,0.0,0.5)
    saveRDS(inis1, file=paste('../results/',date,'/inis.rds', sep = ''))
    wInis = rnorm(N, 0, 1)
    saveRDS(vars, file=paste('../results/',date,'/vars.rds', sep = ''))
  }
  else {
    #load last results (ret object)
    tmp_wInis = readRDS(paste('../results/',date,'/returned.rds', sep=''))
    ncol = length(tmp_wInis$alpha)
    cat(paste("Last round iteration number:", ncol, "\n"))
    wInis = tmp_wInis$w[,ncol]
    inis1 = c(tmp_wInis$alpha[ncol], tmp_wInis$eta[ncol], tmp_wInis$tau2[ncol])
  }
	
  ret1 <- dm_call_cont_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis1, wInis)
  jump_count = get_jump_frequency(ret1, total_iter, N)
  print(jump_count)
  saveRDS(ret1, file = paste('../results/',date,'/returned1.rds', sep = '')) 
  saveRDS(jump_count, file = paste('../results/',date,'/jumps1.rds', sep = ''))

}


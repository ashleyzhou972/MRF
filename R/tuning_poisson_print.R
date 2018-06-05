###Different simulated values and different estimation results


setwd('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/R')
source('../c/dm_call.R')
source('../../functions_for_metropolis.R')
source('../../scaleReductionFactor.R')
source('./simulation_1000.R')
source('./pseudolikelihood.R')
#Total number of iterations
total_iter = 10000
B = 2000
#This observed data are simulated gaussian MRF instead of Poisson
#load observed data
#eta can range from -0.133 to 0.116
###true values:


truevalues = matrix(c(2,-0.1,4,3,0.1,2,1,0.115,7,1,-0.13,2,2,0,5), nrow = 5, byrow = T)
load('./simulated_neighbors_1000_poisson.RData')
net = graph_from_adjacency_matrix(sub_neighbor, mode = "undirected")
LOAD = FALSE
date = "20180605"
#conter = 1
for (counter in 1:5){
  alpha_true= truevalues[counter,1]
  eta_true = truevalues[counter,2]
  tau2_true = truevalues[counter,3]
  filename = paste("./simulated_1000poisson_",alpha_true,"_", eta_true, "_", tau2_true,".RData",sep = "")
  if (LOAD){
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
  
  inis1 = c(5,0.0,0.5)
  ret1 <- dm_call_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis1)
  jump_count = get_jump_frequency(ret1, total_iter, N)
  print(jump_count)
  
  png(filename = paste('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/results/',date,'/sim_1000poisson_plot_', counter, '.png',sep = ""), height = 960, width = 780, pointsize = 14)
  par(mfrow = c(3,1))
  alpha_main = paste("alpha (true= ", alpha_true, ")", sep = "")
  plot_iterations(total_iter, ret1$alpha, inis1[1],bounds_a[1], bounds_a[2],"alpha", alpha_main, alpha_true)
  abline(h = mean(ret1$alpha), col = "#7570b3", lty=2, lwd = 2)
  legend1 = paste("True from simulation:", alpha_true)
  legend3 = paste("Mean MCMC:", round(mean(ret1$alpha),4))
  legend(7000,bounds_a[2]+2,c(legend1, legend3),col = c("#1b9e77","#7570b3"),lty = c(1,2),cex=1.2, xpd=TRUE, lwd = 2)

  eta_main = paste("eta (true= ", eta_true, ")", sep = "")
  plot_iterations(total_iter, ret1$eta, inis1[2],bounds_e[1], bounds_e[2],"eta", eta_main, eta_true)
  abline(h = mean(ret1$eta), col = "#7570b3", lty=2, lwd = 2)
  legend1 = paste("True from simulation:", eta_true)
  legend3 = paste("Mean MCMC:", round(mean(ret1$eta),4))
  legend(7000,bounds_e[2]+0.05,c(legend1, legend3),col = c("#1b9e77","#7570b3"),lty = c(1,2),cex=1.2, xpd=TRUE, lwd = 2)

  tau2_main = paste("tau2 (true= ", tau2_true, ")", sep = "")
  plot_iterations(total_iter, ret1$tau2, inis1[3],bounds_t[1], bounds_t[2],"tau2", tau2_main, tau2_true)
  abline(h = mean(ret1$tau2), col = "#7570b3", lty=2, lwd = 2)
  legend1 = paste("True from simulation:", tau2_true)
  legend3 = paste("Mean MCMC:", round(mean(ret1$tau2),4))
  legend(7000,bounds_t[2]+2,c(legend1, legend3),col = c("#1b9e77","#7570b3"),lty = c(1,2),cex=1.2, xpd=TRUE, lwd = 2)

  dev.off() 
}



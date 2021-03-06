###Different simulated values and different estimation results


setwd('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/R')
source('../c/dm_call.R')
source('../../functions_for_metropolis.R')
#source('../../scaleReductionFactor.R')
source('./simulation_30.R')
source('./pseudolikelihood.R')
library(stats4)
#Total number of iterations
total_iter = 10000
B = 2000
#This observed data are simulated gaussian MRF instead of Poisson
#load observed data
#eta can range from -0.27 to 0.25
###true values:

LOAD = FALSE #if load =1, load simulated y, instead of simulating new y
date = "20180605"
truevalues = matrix(c(2,-0.1,4,3,0.1,2,1,0.2,7,1,-0.2,2,2,0,5), nrow = 5, byrow = T)
load('./simulated_neighbors_30.RData')
for (counter in 1:5){
  if (counter %in% c(1,2,3)) LOAD = TRUE
  else LOAD = FALSE
  alpha_true= truevalues[counter,1]
  eta_true = truevalues[counter,2]
  tau2_true = truevalues[counter,3]
  
  
  filename = paste("./simulated_gaussian_",alpha_true,"_", eta_true, "_", tau2_true,".RData",sep = "")
  if (LOAD){
    load(filename)
    cat("Load simulated y complete\n")
  }else{
    y=simulate_y_gaussian(net,alpha_true,eta_true,tau2_true,10000)
    save(y,file=filename)
    cat("Simulation done and saved\n")
  }
  N = length(y)
  nb_mat_int = as.integer(sub_neighbor)
  #There are four steps of metropolis in each iteration (including double metropolis)
  #When simulating a new value in each of these metropolis steps, variance of the random walk chain can be specified
  vars = c(0.7,0.3,0.5,0.3)
  #first var is ignored in the gaussian case
  #parameters for the prior distributions (uniform), for alpha, eta and tau^2
  bounds_a = c(0,10)
  bounds_e = c(-0.25,0.25)
  bounds_t = c(0,20)
  #inital guess for alpha, beta and tau^2
  
  inis1 = c(5,0.2,0.5)
  ret1 <- dm_call_gaussian_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis1)
  jump_count = get_jump_frequency_gaussian(ret1, total_iter, N)
  print(jump_count)
  

  
  logl_2<-function(alpha, eta, tau2){
    y_arg = y
    neighbor_arg = sub_neighbor
    n = length(y_arg)
    sum_e = 0
    for (i in 1:n){
      mu_i = alpha+eta*neighbor_arg[i,]%*%(y_arg-alpha)
      exp_i = (y_arg[i]-mu_i)^2/tau2
      sum_e = sum_e + exp_i
    }
    l = n*log(tau2)+sum_e
    #This is the negative log-likelihood
    return(l)
  }
  
  #method is BFGS 
  op2<-mle(logl_2,start=list(alpha =1, eta = 1, tau2 = 1))
  print(summary(op2))

  #logLik(op2)
  
  if (counter %in% c(1,2,4)){
    tmax = 20
  } else tmax = 55
  
  tau2_seq = seq(0.5,tmax,0.5)
  tlength = length(tau2_seq)
  ll = rep(NA, tlength)
  alpha_hat = rep(NA, tlength)
  eta_hat = rep(NA, tlength)
  for (tt in tau2_seq){
    results = profile_ll(tt)
    ll[tt*2] = results[[1]]
    op = results[[2]]
    alpha_hat[tt*2] = profile_alpha(op)
    eta_hat[tt*2] = profile_eta(op)
  }
  #Start plotting#
  #date = "20180605"
  
  #three panel MCMC
  png(filename = paste('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/results/',date,'/sim_gaussian_plot_', counter, '.png',sep = ""), height = 960, width = 780, pointsize = 14)
  
  par(mfrow = c(3,1))
  alpha_main = paste("alpha (true= ", alpha_true, ")", sep = "")
  plot_iterations(total_iter, ret1$alpha, inis1[1],bounds_a[1], bounds_a[2],"alpha", alpha_main, alpha_true)
  abline(h = op2@coef[1], col = "#d95f02", lty=2, lwd = 2)
  abline(h = mean(ret1$alpha), col = "#7570b3", lty=2, lwd = 2)
  legend1 = paste("True from simulation:", alpha_true)
  legend2 = paste("M_Pseudo_LE:", round(op2@coef[1],4))
  legend3 = paste("Mean MCMC:", round(mean(ret1$alpha),4))
  legend(7000,bounds_a[2]+2,c(legend1, legend2, legend3),col = c("#1b9e77","#d95f02","#7570b3"),lty = c(1,2,2),cex=1.2, xpd=TRUE, lwd = 2)
  
  eta_main = paste("eta (true= ", eta_true, ")", sep = "")
  plot_iterations(total_iter, ret1$eta, inis1[2],bounds_e[1], bounds_e[2],"eta", eta_main, eta_true)
  abline(h = op2@coef[2], col = "#d95f02", lty=2, lwd = 2)
  abline(h = mean(ret1$eta), col = "#7570b3", lty=2, lwd = 2)
  legend1 = paste("True from simulation:", eta_true)
  legend2 = paste("M_Pseudo_LE:", round(op2@coef[2],4))
  legend3 = paste("Mean MCMC:", round(mean(ret1$eta),4))
  legend(7000,bounds_e[2]+0.12,c(legend1, legend2, legend3),col = c("#1b9e77","#d95f02","#7570b3"),lty = c(1,2,2),cex=1.2, xpd=TRUE, lwd = 2)
  
  tau2_main = paste("tau2 (true= ", tau2_true, ")", sep = "")
  plot_iterations(total_iter, ret1$tau2, inis1[3],bounds_t[1], bounds_t[2],"tau2", tau2_main, tau2_true)
  abline(h = op2@coef[3], col = "#d95f02", lty=2, lwd = 2)
  abline(h = mean(ret1$tau2), col = "#7570b3", lty=2, lwd = 2)
  legend1 = paste("True from simulation:", tau2_true)
  legend2 = paste("M_Pseudo_LE:", round(op2@coef[3],4))
  legend3 = paste("Mean MCMC:", round(mean(ret1$tau2),4))
  legend(7000,bounds_t[2]+2,c(legend1, legend2, legend3),col = c("#1b9e77","#d95f02","#7570b3"),lty = c(1,2,2),cex=1.2, xpd=TRUE, lwd = 2)
  
  dev.off()
  
  #single panel profile likelihood
  png(filename = paste('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/results/',date,'/sim_gaussian_profile_single', counter, '.png',sep = ""), height = 400, width = 780, pointsize = 14)
  
  par()
  plot(x=tau2_seq,y = ll, xlab = "tau2", ylab = "log likelihood", main = "profile likelihood on tau2", pch = 19, pwd = 1.2,cex.lab  = 1, cex.main = 1)
  abline(v = truevalues[counter,3], col = 3, lwd = 2)
  abline(v = tau2_seq[which(ll==max(ll))], col = 4, lty = 2, lwd = 2)
  legend1 = paste("true simulated tau2:", truevalues[counter,3])
  legend2 = paste("tau2 estimate from pseudolikelihood:", tau2_seq[which(ll==max(ll))])
  legend("bottomright", legend = c(legend1,legend2),col = c(3,4), lty = c(1,2),lwd =2, cex = 1)
  dev.off()
  
  #three panel profile likelihood
  png(filename = paste('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/results/',date,'/sim_gaussian_profile_', counter, '.png',sep = ""), height = 1024, width = 780, pointsize = 14)
  
  par(mfrow = c(3,1), xpd = T)
  plot(x=tau2_seq,y = ll, xlab = "tau2", ylab = "log likelihood", main = "profile likelihood on tau2", cex.lab  = 1.7, cex.main = 1.5)
  abline(v = truevalues[counter,3], col = 3, lwd = 2)
  abline(v = tau2_seq[which(ll==max(ll))], col = 4, lty = 2, lwd = 2)
  legend1 = paste("true simulated tau2:", truevalues[counter,3])
  legend2 = paste("tau2 estimate from pseudolikelihood:", tau2_seq[which(ll==max(ll))])
  legend("bottomright", legend = c(legend1,legend2),col = c(3,4), lty = c(1,2),lwd =2, cex = 1.7)
  
  plot(x=tau2_seq,y = alpha_hat, xlab = "tau2", ylab = "alpha_hat", cex.lab  = 1.7, cex.main = 1.5)
  abline(v = truevalues[counter,3], col = 3, lwd = 2)
  abline(v = tau2_seq[which(ll==max(ll))], col = 4, lty = 2, lwd = 2)
  plot(x=tau2_seq,y = eta_hat , xlab = "tau2", ylab = "eta_hat", cex.lab  = 1.7, cex.main = 1.5)
  abline(v = truevalues[counter,3], col = 3, lwd = 2)
  abline(v = tau2_seq[which(ll==max(ll))], col = 4, lty = 2, lwd = 2)
  dev.off()
}


###Different simulated values and different estimation results


setwd('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/R')
source('../c/dm_call.R')
source('../../functions_for_metropolis.R')
source('../../scaleReductionFactor.R')
source('./simulation_30.R')
#Total number of iterations
total_iter = 20000
B = 2000
#This observed data are simulated gaussian MRF instead of Poisson
#load observed data
#eta can range from -0.27 to 0.25
###true values:


truevalues = matrix(c(2,-0.1,4,3,0.1,2,1,0.2,7,1,-0.2,2,2,0,5), nrow = 5, byrow = T)

for (counter in 1:5){
  alpha_true= truevalues[counter,1]
  eta_true = truevalues[counter,2]
  tau2_true = truevalues[counter,3]
  y=simulate_y_gaussian(net,alpha_true,eta_true,tau2_true,20000)
  #save(y,file='./simulated_y_gaussian.RData')
  #load('./simulated_y_gaussian.RData')
  load('./simulated_neighbors_30.RData')
  
  N = length(y)
  nb_mat_int = as.integer(sub_neighbor)
  #There are four steps of metropolis in each iteration (including double metropolis)
  #When simulating a new value in each of these metropolis steps, variance of the random walk chain can be specified
  vars = c(0.7,0.3,0.5,0.3)
  #first var is ignored in the gaussian case
  #parameters for the prior distributions (uniform), for alpha, eta and tau^2
  bounds_a = c(0,10)
  bounds_e = c(-0.27,0.25)
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
    sum = 0
    for (i in 1:n){
      mu_i = alpha+eta*neighbor_arg[i,]%*%(y_arg-alpha)
      exp_i = (y_arg[i]-mu_i)^2
      sum = sum + exp_i
    }
    l = n*log(tau2)+sum/tau2
    #This is the negative log-likelihood
    return(l)
  }
  
  #method is BFGS 
  op2<-mle(logl_2,start=list(alpha =1, eta = 1, tau2 = 1))
  summary(op2)
  logLik(op2)
  
  png(filename = paste('./sim_all_plot_', counter, '.png',sep = ""), height = 960, width = 780, pointsize = 14)
  
  par(mfrow = c(3,1))
  alpha_main = paste("alpha (true= ", alpha_true, ")", sep = "")
  plot_iterations(total_iter, ret1$alpha, inis1[1],bounds_a[1], bounds_a[2],"alpha", alpha_main, alpha_true)
  abline(h = op2@coef[1], col = "#d95f02", lty=2, lwd = 2)
  abline(h = mean(ret1$alpha), col = "#7570b3", lty=2, lwd = 2)
  legend1 = paste("True from simulation:", alpha_true)
  legend2 = paste("M_Pseudo_LE:", round(op2@coef[1],4))
  legend3 = paste("Mean MCMC:", round(mean(ret1$alpha),4))
  legend(12000,bounds_a[2]+2,c(legend1, legend2, legend3),col = c("#1b9e77","#d95f02","#7570b3"),lty = c(1,2,2),cex=1.2, xpd=TRUE, lwd = 2)
  
  eta_main = paste("eta (true= ", eta_true, ")", sep = "")
  plot_iterations(total_iter, ret1$eta, inis1[2],bounds_e[1], bounds_e[2],"eta", eta_main, eta_true)
  abline(h = op2@coef[2], col = "#d95f02", lty=2, lwd = 2)
  abline(h = mean(ret1$eta), col = "#7570b3", lty=2, lwd = 2)
  legend1 = paste("True from simulation:", eta_true)
  legend2 = paste("M_Pseudo_LE:", round(op2@coef[2],4))
  legend3 = paste("Mean MCMC:", round(mean(ret1$eta),4))
  legend(12000,bounds_e[2]+0.12,c(legend1, legend2, legend3),col = c("#1b9e77","#d95f02","#7570b3"),lty = c(1,2,2),cex=1.2, xpd=TRUE, lwd = 2)
  
  tau2_main = paste("tau2 (true= ", tau2_true, ")", sep = "")
  plot_iterations(total_iter, ret1$tau2, inis1[3],bounds_t[1], bounds_t[2],"tau2", tau2_main, tau2_true)
  abline(h = op2@coef[3], col = "#d95f02", lty=2, lwd = 2)
  abline(h = mean(ret1$tau2), col = "#7570b3", lty=2, lwd = 2)
  legend1 = paste("True from simulation:", tau2_true)
  legend2 = paste("M_Pseudo_LE:", round(op2@coef[3],4))
  legend3 = paste("Mean MCMC:", round(mean(ret1$tau2),4))
  legend(12000,bounds_t[2]+2,c(legend1, legend2, legend3),col = c("#1b9e77","#d95f02","#7570b3"),lty = c(1,2,2),cex=1.2, xpd=TRUE, lwd = 2)
  
  dev.off()
}



logl<-function(params, y, neighbor){
  alpha = params[1]
  eta = params[2]
  tau2 = params[3]
  n = length(y)
  sum = 0
  for (i in 1:n){
    mu_i = alpha+eta*neighbor[i,]%*%(y-alpha)
    exp_i = (y[i]-mu_i)^2
    sum = sum + exp_i
  }
  l = n*log(tau2)+sum/tau2
  #This is the negative log-likelihood
  return(l)
}



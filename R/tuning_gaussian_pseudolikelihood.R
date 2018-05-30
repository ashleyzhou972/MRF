setwd('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/R')
library(stats4)
source('./simulation_30.R')

profile_ll<-function(tau2){
  op2<-mle(logl_2,start=list(alpha =1, eta = 1),fixed = list(tau2 = tau2))
  #summary(op2)
  r = logLik(op2)
  return(r)
}

profile_alpha<-function(tau2){
  op2<-mle(logl_2,start=list(alpha =1, eta = 1),fixed = list(tau2 = tau2))
  #summary(op2)
  r = op2@coef[1]
  return(r)
}


profile_eta<-function(tau2){
  op2<-mle(logl_2,start=list(alpha =1, eta = 1),fixed = list(tau2 = tau2))
  #summary(op2)
  r = op2@coef[2]
  return(r)
}


#This observed data are simulated gaussian MRF instead of Poisson
#load observed data
#eta can range from -0.25 to 0.25
###true values:


truevalues = matrix(c(2,-0.1,4,3,0.1,2,1,0.2,7,1,-0.2,2,2,0,5), nrow = 5, byrow = T)
load('./simulated_neighbors_30.RData')
for (counter in 1:5){
  alpha_true= truevalues[counter,1]
  eta_true = truevalues[counter,2]
  tau2_true = truevalues[counter,3]
  filename = paste("./simulated_gaussian_",alpha_true,"_", eta_true, "_", tau2_true,".RData",sep = "")
  load(filename)
  cat("Load simulated y complete\n")
  #save(y,file='./simulated_y_gaussian.RData')
  #load('./simulated_y_gaussian.RData')
  
  logl_2<-function(alpha, eta, tau2){
    ##This is function to calculate negative loglikelihood
    y_arg = y
    neighbor_arg = sub_neighbor
    n = length(y_arg)
    sum = 0
    for (i in 1:n){
      mu_i = alpha+eta*neighbor_arg[i,]%*%(y_arg-alpha)
      exp_i = (y_arg[i]-mu_i)^2/tau2
      sum = sum + exp_i
    }
    l = n*log(tau2)+sum
    #This is the negative log-likelihood
    return(l)
  }
  
  if (counter %in% c(1,2,4)){
    tmax = 20
  }
  else tmax = 55
  tau2_seq = seq(0.5,tmax,0.5)
  tlength = length(tau2_seq)
  ll = rep(NA, tlength)
  alpha_hat = rep(NA, tlength)
  eta_hat = rep(NA, tlength)
  for (t in tau2_seq){
    ll[t*2] = profile_ll(t)
    alpha_hat[t*2] = profile_alpha(t)
    eta_hat[t*2] = profile_eta(t)
  }
  
  #three panel
  png(filename = paste('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/results/20180529','/sim_gaussian_profile_', counter, '.png',sep = ""), height = 1024, width = 780, pointsize = 14)
  
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
  
  #single panel
  png(filename = paste('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/results/20180529','/sim_gaussian_profile_single', counter, '.png',sep = ""), height = 400, width = 780, pointsize = 14)
  
  par()
  plot(x=tau2_seq,y = ll, xlab = "tau2", ylab = "log likelihood", main = "profile likelihood on tau2", pch = 19, pwd = 1.2,cex.lab  = 1, cex.main = 1)
  abline(v = truevalues[counter,3], col = 3, lwd = 2)
  abline(v = tau2_seq[which(ll==max(ll))], col = 4, lty = 2, lwd = 2)
  legend1 = paste("true simulated tau2:", truevalues[counter,3])
  legend2 = paste("tau2 estimate from pseudolikelihood:", tau2_seq[which(ll==max(ll))])
  legend("bottomright", legend = c(legend1,legend2),col = c(3,4), lty = c(1,2),lwd =2, cex = 1)
  dev.off()
}


#################################################################

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

logl_2<-function(alpha, eta, tau2){
  ##This is function to calculate negative loglikelihood
  y_arg = y
  neighbor_arg = sub_neighbor
  n = length(y_arg)
  sum = 0
  for (i in 1:n){
    mu_i = alpha+eta*neighbor_arg[i,]%*%(y_arg-alpha)
    exp_i = (y_arg[i]-mu_i)^2/tau2
    sum = sum + exp_i
  }
  l = n*log(tau2)+sum
  #This is the negative log-likelihood
  return(l)
}


#method is Nelder-Mead
# op1<-optim(truevalues[counter,], logl,y=y, neighbor = sub_neighbor)
# #method is BFGS 
# op2<-mle(logl_2,start=list(alpha =1, eta = -0.1, tau2 = 1))
# summary(op2)
# logLik(op2)

###########################fix tau2################
logl_3<-function(params, tau2, y, neighbor){
  alpha = params[1]
  eta = params[2]
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
# 
# op3<-optim(truevalues[counter,c(1,2)],logl_3, tau2 = 4, y = y, neighbor = sub_neighbor)
# op3$value
# op4<-mle(logl_2, start=list(alpha=1, eta=-0.1),fixed = list(tau2 = 4))
# summary(op4)
# logLik(op4)




#####################################################################
#profile likelihood on tau2
#####################################################################






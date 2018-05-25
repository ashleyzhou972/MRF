library(stats4)
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
#method is Nelder-Mead
op1<-optim(truevalues[counter,], logl,y=y, neighbor = sub_neighbor)
#method is BFGS 
op2<-mle(logl_2,start=list(alpha =1, eta = -0.1, tau2 = 1))
summary(op2)
logLik(op2)

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

op3<-optim(truevalues[counter,c(1,2)],logl_3, tau2 = 4, y = y, neighbor = sub_neighbor)
op3$value
op4<-mle(logl_2, start=list(alpha=1, eta=-0.1),fixed = list(tau2 = 4))
summary(op4)
logLik(op4)

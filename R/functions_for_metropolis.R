library(MASS)
#library(mvtnorm)
#library(Matrix)
library(igraph, quietly = T,warn.conflicts = F)

find_con_cliques<-function(neighbor_matrix,k){
  #There should be three outcomes to this
  #k=1: Find concliques, and each conclique IDs are in a list
  #k=2: Do not find concliques, and return each individual observation as a conclique
  #k=3: Do not find concliques, and return the observations as a whole (marginal joint)
  #k=3 has a different target distribution for Metropollis than the first two
  #Use k=1 only with simulated (small) dataset
  N = dim(neighbor_matrix)[1]
  if (k==1){
    comp_matrix = matrix(1,nrow=100,ncol=100)-neighbor_matrix
    c_network = graph_from_adjacency_matrix(comp_matrix,mode="undirected",diag=F)
    quantile(degree(c_network))
    cliques =names(largest_cliques(c_network)[[1]])
    #The largest clique has 52 members
    
    clique1_id = which(rownames(neighbor_matrix)%in%cliques)
    clique1_names = rownames(neighbor_matrix)[clique1_id]
    #check
    #sum(cliques %in% clique1_names)
    clique2_names = rownames(neighbor_matrix)[-clique1_id]
    clique2_id = which(rownames(neighbor_matrix)%in%clique2_names)
    return(list(conclique1=clique1_id,conclique2=clique2_id))
  }else if (k==2){
    return(as.list(c(1:N)))
  }else if (k==3){
    return(list(c(1:N)))
  }
  
}
data_density<-function(w,y){
  #f(y_i|w_i) = \dfrac{1}{y_i!}exp(-e^{w_i{}+w_iy_i+w_i)
  #y is the observed data
  #both in vector
  #should be log density
  f = -sum(lfactorial(y))-sum(exp(w))+crossprod(w,y)
  return(f)
}

normal_density<-function(y,mean,covariance){
  if (length(y)==1){
    d = dnorm(y,mean,sd=covariance,log=TRUE)
  }else{
    #This is the multivariate version
    d=dmvnorm(y,mean,sigma=covariance,log=TRUE)  
  }
  return(d)
}
proposal<-function(old,new,prior_function,data_function,q_function){
  #both parameters (old and new) should be vectors
  #The last three parameters should be functions with one or two parameters
  #All densities should be in log form!!!!! (20171003)
  nominator = prior_function(new)+data_function(new)+q_function(old,new)
  denominator = prior_function(old)+data_function(old)+q_function(new,old)
  lalpha = nominator-denominator
  #cat(lalpha,'\n')
  ######################################################################3
  #TODO
  #What if the denominator is 0?!
  if (is.na(lalpha)) lalpha=-Inf
  alpha=exp(lalpha)
  if (alpha<=1) return(alpha)
  else return(1)
}
sigma_calc<-function(neighbor_matrix,eta,tau){
  n = dim(neighbor_matrix)[1]
  i_c = diag(n)-eta*neighbor_matrix
  i_c_inv = solve(i_c)
  sig = crossprod(i_c_inv,diag(n)*tau)
  #(I-C)^{-1} is a symmetric matrix
  #so no need to transpose
  #Diagonal should be positive
  return(sig)
}


##In all current cases q_functions are the same for both sides of the fraction
##So we define a generic q_function to be plugged in for the function "proposal"
q_function_generic<-function(old,new){
  return(1)
}

scaleReductionFactor<-function(param){
  #This is for ONE parameter in the model
  #param should be a matrix of n rows and M columns
  #n is the number of iterations so far
  #M is the number of starting values (chains) we tried
  n = nrow(param)
  #M = ncol(param)
  sample_var = apply(param,MARGIN=2,FUN=var)
  within_seq_var = mean(sample_var)
  sample_mean = colMeans(param)
  between_seq_var = n*var(sample_mean)
  var_hat = (n-1)/n*within_seq_var+1/n*between_seq_var
  r = sqrt(var_hat/within_seq_var)
  return(r)
}


plot_iterations_srf<-function(Ti,B,alpha,alpha0,ylim1,ylim2,name){
  plot(x=c(1:Ti),y=alpha[1:Ti,1],type='l',xlab = 'iteration',ylab=name,ylim = c(ylim1,ylim2))
  lines(x=c(1:Ti),y=alpha[1:Ti,2],type='l',col=2)
  lines(x=c(1:Ti),y=alpha[1:Ti,3],type='l',col=3)
  lines(x=c(1:Ti),y=alpha[1:Ti,4],type='l',col=4)
  lines(x=c(1:Ti),y=alpha[1:Ti,5],type='l',col=5)
  abline(v=B,col= "#7fc97f")
  legend("topright",legend=paste('Starting value:', alpha0),col=c(1,2,3,4,5),lty=1)
}

plot_iterations<-function(Ti,B,alpha,alpha0,ylim1,ylim2,name,main, true){
  plot(x=c(1:Ti),y=alpha[1:Ti],type='l',xlab = 'iteration',ylab=name,ylim = c(ylim1,ylim2), main = main, bty = "L", lwd = 2)
  #abline(v=B,col=6)
  abline(h = true, col = "#1b9e77")
  #legend("topright",legend=paste('Starting value:', alpha0),col=c(1,2,3,4,5),lty=1)
}
get_posterior_distribution_srf<-function(Ti, B,alpha){
  return(c(alpha[B:Ti,]))
}

get_posterior_distribution<-function(Ti, B,alpha){
  return(alpha[B:Ti])
}

assemble_rets<-function(ret1, ret2) {
#This function concatenates returned objects from two results
#ret2 should have initial values as the last iteration of ret1
  newret = list()
  newret$w = cbind(ret1$w, ret2$w)
  newret$alpha = c(ret1$alpha, ret2$alpha)
  newret$eta = c(ret1$eta, ret2$eta)
  newret$tau2 = c(ret1$tau2, ret2$tau2)
  return(newret)
}

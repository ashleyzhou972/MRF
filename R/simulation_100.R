#!/usr/bin/env Rscript
#Simulate a smaller set of data
#n=100
#neighborhood informatino randomly drawn from conclique (belongs to two concliques)

#Since the dependence parameter eta should be bounded by characteristics of the adjacency matrix
#The simulated data should have a valid eta

#library(MASS)
#library(mvtnorm)
#library(RcppEigen)
library(Matrix)
#source('./read.R')
library(igraph)

#Updated in 20180122:
#because trace is zero (sum of diagonal values), the max and min eigenvalues are opposite signs
#so parameter space of eta is h_1^{-1}<eta<h_n^{-1}, where h_1<h_2<...<h_n are eigenvalues

#The larger the absolute eigenvalues, the smaller the sample space
#How does that connect to eigencentrality?

#Test for different network characteristics
#Simulate a network with 100 nodes
#To preserve the network characteristics of the original model
#scale free network

#characteristics of the original gene network
#net = graph_from_adjacency_matrix(c,mode="undirected")
#dd = degree(net)
#hist(dd,prob=T,main="Histogram of Degrees of HiC Network",xlab="Degree")
#summary(dd)
##node to edge ratio
#vcount(net)/ecount(net)
##power-law (scale free)
#pl<-fit_power_law(dd)
###alpha=4
###p=0.75


#simulate smaller network with the same characteristics
##type1, erdos-renyi
##same v/e ratio
##BELOW IS JUST EXAMPLE
#rg1<-erdos.renyi.game(100,384,type='gnm')
##type, power-law with fitness
#rg2<-sample_fitness_pl(100,384,4)
#hist(degree(rg2),prob=T)
#fit_power_law(degree(rg2))
##type3, power-law with preferential attachment (Barabasi-Albert)
#rg3<-sample_pa(100,m=2,zero.appeal = 2,directed=F)
#hist(degree(rg3),prob=T)
#fit_power_law(degree(rg3))

generate_random_graph<-function(type,v,e,power){
  if (type==1){
    net = erdos.renyi.game(v,e,type="gnm")
  }
  else if (type==2){
    net = sample_fitness_pl(v,e,power)
  }
  else if (type==3){
    net = sample_pa(v,m=2,zero.appeal = 2,directed=F)
  }
  return(net)
}

get_eigen_interval<-function(net){
  adj_mat = as_adj(net,type="both",sparse=F)
  evalues = eigen(adj_mat)$values
  e_max = max(evalues)
  e_min = min(evalues)
  #e_max and e_min should be opposite signs
  r = 1/e_max-1/e_min
  return(list(1/e_min,1/e_max,r))
}

#eigen values for the original matrix
#max=44.1972
#min=-12.6505
#eigen_interval = (-0.07905,0.02263)=0.1016


simulate_y<-function(net,alpha,eta,tau2,M){
  #Simulate w's
  #using Gibbs
  n = vcount(net)
  sub_neighbor = as_adj(net,type="both",sparse=F)
  for (t in 1:M){
    w = rep(0,n)
    for (i in 1:n){
      mu = alpha+eta*sub_neighbor[i,]%*%(w-alpha)
      w[i] = rnorm(1,mu,tau2)
    }
  }
  lambda = exp(w)
  y = rep(NA,n)
  for (i in 1:n){
    y[i] = rpois(1,lambda[i])
  }
  return(y)
}


#####main######
# 
# set.seed(2)
# net = generate_random_graph(3,100,384,4)
# hist(degree(net))
# #power-law distribution
# components(net)
# #net is connected 
# get_eigen_interval(net)
# #eta can range from -0.199 to 0.154
# 
# ###true values:
# #alpha=2
# #eta = 0.1
# #tau2 = 2
# y=simulate_y(net,2,0.1,2,10000)
# 
# 
# save(y,file='./simulated_y__new.RData')
# sub_neighbor = as_adj(net,type="both",sparse=F)
# save(sub_neighbor,file = './simulated_neighbors_new.RData')
# 
# 





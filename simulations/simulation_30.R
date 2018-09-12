#!/usr/bin/env Rscript
#Simulate a smaller set of data
#n=30
#The data simulated will be gaussian

#Since the dependence parameter eta should be bounded by characteristics of the adjacency matrix
#The simulated data should have a valid eta

#setwd('/home/nzhou/hic/IMR90/work')
#library(MASS)
#library(mvtnorm)
#library(RcppEigen)
library(Matrix)
#source('./read.R')
library(igraph,warn.conflicts = F,quietly = T)

generate_four_nearest_neighbor_matrix<-function(v,k){
  #v is number of nodes (be an even number please)
  #k is how v is divided into a rectangular field
  #for example if v = 30
  #k = 5 (6*5 = 30)
  #k is length of column
  #1  6 11 16 21 26#
  #2  7 12 17 22 27#
  #3  8 13 18 23 28#
  #4  9 14 19 24 29#
  #5 10 15 20 25 30#
  mat = matrix(0, nrow = v, ncol = v)
  #four vertex cases (1, k, 1+v-k and v)
  #four edge cases ([1,k], [v-k+1,v],[i%%k==1], [i%%k==1])
  #rest are internal cases
  pointy = c(1,k,1+v-k,v)
  edgeL = c(1:k)
  edgeR = c((v-k+1):v)
  rlength = v/k
  edgeT = rep(NA, rlength)
  edgeB = rep(NA, rlength)
  for (j in 1:rlength){
    edgeT[j] = 1+(j-1)*k
    edgeB[j] = j*k
  }
  edges = c(edgeL, edgeR, edgeT, edgeB)
  for (i in 1:v){
    if (i %in% edges){
      if (i %in% edgeL){
        mat[i,i+k] = 1
        mat[i,i-k+v] = 1
        if (!(i %in% pointy)){
          mat[i,i+1] = 1
          mat[i,i-1] = 1
        }
      }
      if (i %in% edgeR){
        mat[i,i+k-v] = 1
        mat[i,i-k] = 1
        if (!(i %in% pointy)){
          mat[i,i+1] = 1
          mat[i,i-1] = 1
        }
      }
      if (i %in% edgeT){
        mat[i,i+1] = 1
        mat[i,i-1+k] = 1
        if (!(i %in% pointy)){
          mat[i,i+k] = 1
          mat[i,i-k] = 1
        }
      }
      if (i %in% edgeB){
        mat[i,i+1-k] = 1
        mat[i,i-1] = 1
        if (!(i %in% pointy)){
          mat[i,i+k] = 1
          mat[i,i-k] = 1
        }
      }
    }
    else {
      mat[i,i+1] = 1
      mat[i,i-1] = 1
      mat[i,i+k] = 1
      mat[i,i-k] = 1
    }
  }
  return(mat)
}




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


simulate_y_poisson<-function(net,alpha,eta,tau2,M){
  #using Gibbs
  n = vcount(net)
  sub_neighbor = as_adj(net,type="both",sparse=F)
  w = rep(0,n)
  for (t in 1:M){
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


simulate_y_gaussian<-function(net,alpha,eta,tau2,M){
  #using Gibbs
  n = vcount(net)
  sub_neighbor = as_adj(net,type="both",sparse=F)
  #above has no diagonal values
  #sub_neighbor = sub_neighbor + diag(n)
  w = rep(0,n)
  for (t in 1:M){
    cat('iteration',t,'\n')
    for (i in 1:n){
      mu = alpha+eta*sub_neighbor[i,]%*%(w-alpha)
      w[i] = rnorm(1,mu,tau2)
      #print(w)
    }
  }
  return(w)
}

#####main######

set.seed(2)
net = graph_from_adjacency_matrix(generate_four_nearest_neighbor_matrix(900,30),mode = "undirected")
#quantile(degree(net))
#power-law distribution
#components(net)
#net is connected 
print(get_eigen_interval(net)[1:2])
#eta can range from -0.27 to 0.25



# 
# save(y,file='./simulated_y_gaussian.RData')
# sub_neighbor = as_adj(net,type="both",sparse=F)
# save(sub_neighbor,file = './simulated_neighbors_30.RData')
#






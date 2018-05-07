############################################
# R wrapper for double metropolis
############################################

setwd('/home/nzhou/hic/IMR90/work/dm/MRF_HIC_GE/c')
dm_call_wrapper<-function(total_iter,y, nb_mat,vars, bounds_a, bounds_e, bounds_t, inis){
	if (!is.numeric(total_iter) || !is.numeric(y) || !is.numeric(nb_mat) || !is.numeric(vars) || !is.numeric(bounds_a) || !is.numeric(bounds_e) || !is.numeric(bounds_t) || !is.numeric(inis)){
		stop("input data not numeric\n")
	}
	if (!is.loaded("double_metropolis", PACKAGE="dm_call")){
		print("loading")
		dyn.load("dm_call.so")
	}
	#print(is.loaded("double_metropolis", PACKAGE="dm_call"))
	ret = .Call("double_metropolis", T_in = as.integer(total_iter), y_in = as.double(y),neighbor_in = as.integer(nb_mat),vars_in =as.double(vars), bounds_alpha = as.double(bounds_a), bounds_eta = as.double(bounds_e), bounds_tau2 = as.double(bounds_t),as.double(inis))
}

#Total number of iterations
total_iter = 2000
#load observed data
load('./simulated_y__new.RData')
load('./simulated_neighbors_new.RData')
nb_mat_int = as.integer(sub_neighbor)
#There are four steps of metropolis in each iteration (including double metropolis)
#When simulating a new value in each of these metropolis steps, variance of the random walk chain can be specified
vars = c(0.1,0.2,0.3,0.4)
#parameters for the prior distributions (uniform), for alpha, eta and tau^2
bounds_a = c(0,5)
bounds_e = c(-0.15,0.15)
bounds_t = c(0,5)
#inital guess for alpha, beta and tau^2
inis = c(2,0.1,2)
ret <- dm_call_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis)
#results
alpha = ret[[2]]
eta = ret[[3]]
tau2 = ret[[4]]


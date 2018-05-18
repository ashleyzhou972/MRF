############################################
# R wrapper for double metropolis
############################################

dm_call_wrapper<-function(total_iter,y, nb_mat,vars, bounds_a, bounds_e, bounds_t, inis){
	if (!is.numeric(total_iter) || !is.numeric(y) || !is.numeric(nb_mat) || !is.numeric(vars) || !is.numeric(bounds_a) || !is.numeric(bounds_e) || !is.numeric(bounds_t) || !is.numeric(inis)){
		stop("input data not numeric\n")
	}
	if (!is.loaded("double_metropolis", PACKAGE="dm_call")){
		print("loading")
		dyn.load("../c/dm_call.so")
	}
	#print(is.loaded("double_metropolis", PACKAGE="dm_call"))
	ret = .Call("double_metropolis", T_in = as.integer(total_iter), y_in = as.double(y),neighbor_in = as.integer(nb_mat),vars_in =as.double(vars), bounds_alpha = as.double(bounds_a), bounds_eta = as.double(bounds_e), bounds_tau2 = as.double(bounds_t),as.double(inis))
	return(list(w = ret[[1]], alpha = ret[[2]], eta = ret[[3]], tau2 = ret[[4]], jump_count = ret[[5]]))
}



get_jump_frequency<-function(ret, Ti, N){
  #ret is a list of returned values from c
  #Ti is total number of iterations
  #N is data size
  
  jumpcount = ret[[5]]
  rate_w = jumpcount[1]/(N*Ti)
  rate_alpha = jumpcount[2]/Ti
  rate_eta = jumpcount[3]/Ti
  rate_tau2 = jumpcount[4]/Ti
  
  return(list(w = rate_w, alpha = rate_alpha, eta = rate_eta, tau2 = rate_tau2))
}
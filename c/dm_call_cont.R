############################################
# R wrapper for double metropolis
############################################

dm_call_wrapper<-function(total_iter,y, nb_mat,vars, bounds_a, bounds_e, bounds_t, inis){
	if (!is.numeric(total_iter) || !is.numeric(y) || !is.numeric(nb_mat) || !is.numeric(vars) || !is.numeric(bounds_a) || !is.numeric(bounds_e) || !is.numeric(bounds_t) || !is.numeric(inis)){
		stop("input data not numeric\n")
	}
	else {
	  cat("loading package...\n")
	  dyn.load("../c/dm_call.so")
	  cat("loaded status:", is.loaded("double_metropolis", PACKAGE="dm_call"),"\n")
	}
	#print(is.loaded("double_metropolis", PACKAGE="dm_call"))
	ret = .Call("double_metropolis", T_in = as.integer(total_iter), y_in = as.double(y),neighbor_in = as.integer(nb_mat),vars_in =as.double(vars), bounds_alpha = as.double(bounds_a), bounds_eta = as.double(bounds_e), bounds_tau2 = as.double(bounds_t),as.double(inis))
	return(list(w = ret[[1]], alpha = ret[[2]], eta = ret[[3]], tau2 = ret[[4]], jump_count = ret[[5]]))
}


dm_call_cont_wrapper<-function(total_iter,y, nb_mat,vars, bounds_a, bounds_e, bounds_t, inis, inis_w){
	if (!is.numeric(total_iter) || !is.numeric(y) || !is.numeric(nb_mat) || !is.numeric(vars) || !is.numeric(bounds_a) || !is.numeric(bounds_e) || !is.numeric(bounds_t) || !is.numeric(inis) || !is.numeric(inis_w)){
		stop("input data not numeric\n")
	}
	else {
	  cat("loading package...\n")
	  dyn.load("../c/dm_call_cont.so")
	  cat("loaded status:", is.loaded("double_metropolis_cont", PACKAGE="dm_call_cont"),"\n")
	}
	#print(is.loaded("double_metropolis", PACKAGE="dm_call"))
	ret = .Call("double_metropolis_cont", T_in = as.integer(total_iter), y_in = as.double(y),neighbor_in = as.integer(nb_mat),vars_in =as.double(vars), bounds_alpha = as.double(bounds_a), bounds_eta = as.double(bounds_e), bounds_tau2 = as.double(bounds_t), inis=as.double(inis), wInitial=as.double(inis_w))
	return(list(w = ret[[1]], alpha = ret[[2]], eta = ret[[3]], tau2 = ret[[4]], jump_count = ret[[5]]))
}


dm_call_gaussian_wrapper<-function(total_iter,y, nb_mat,vars, bounds_a, bounds_e, bounds_t, inis){
	if (!is.numeric(total_iter) || !is.numeric(y) || !is.numeric(nb_mat) || !is.numeric(vars) || !is.numeric(bounds_a) || !is.numeric(bounds_e) || !is.numeric(bounds_t) || !is.numeric(inis)){
		stop("input data not numeric\n")
	}
	else{
		cat("loading package...\n")
		dyn.load("../c/dm_call_gaussian.so")
		cat("loaded status:", is.loaded("double_metropolis_gaussian", PACKAGE="dm_call_gaussian"),"\n")
	}
	ret = .Call("double_metropolis_gaussian", T_in = as.integer(total_iter), y_in = as.double(y),neighbor_in = as.integer(nb_mat),vars_in =as.double(vars), bounds_alpha = as.double(bounds_a), bounds_eta = as.double(bounds_e), bounds_tau2 = as.double(bounds_t),as.double(inis))
	return(list(alpha = ret[[1]], eta = ret[[2]], tau2 = ret[[3]], jump_count = ret[[4]]))
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



get_jump_frequency_gaussian<-function(ret, Ti, N){
  #ret is a list of returned values from c
  #Ti is total number of iterations
  #N is data size
  
  jumpcount = ret[[4]]
  rate_alpha = jumpcount[1]/Ti
  rate_eta = jumpcount[2]/Ti
  rate_tau2 = jumpcount[3]/Ti
  
  return(list(alpha = rate_alpha, eta = rate_eta, tau2 = rate_tau2))
}

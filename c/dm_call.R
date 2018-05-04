############################################
#R wrapper for double metropolis
############################################

setwd('/home/nzhou/stat580/finalProject/c/')
dm_call_wrapper<-function(total_iter,y, nb_mat,vars, bounds_a, bounds_e, bounds_t){
	if (!is.numeric(total_iter) || !is.numeric(y) || !is.numeric(nb_mat) || !is.numeric(vars) || !is.numeric(bounds_a) || !is.numeric(bounds_e) || !is.numeric(bounds_t)){
		stop("input data not numeric\n")
	}
	if (!is.loaded("double_metropolis", PACKAGE="dm_call")){
		print("loading")
		dyn.load("dm_call.so")
	}
	#print(is.loaded("double_metropolis", PACKAGE="dm_call"))
	ret = .Call("double_metropolis", T_in = as.integer(total_iter), y_in = as.integer(y),neighbor_in = as.integer(nb_mat),vars_in =as.double(vars), bounds_alpha = as.double(bounds_a), bounds_eta = as.double(bounds_e), bounds_tau2 = as.double(bounds_t))
}

total_iter = 100
y = c(1,2,3,4,5,6,7,8,9)
set.seed(2)
nb_mat = matrix(rbinom(72,1,.5),nrow = 9)
nb_mat1 = cbind(c(1,0,0,0,0,0,0,0,1),nb_mat)
nb_mat_int = as.integer(nb_mat1)
#is.integer(nb_mat_int)
#is.vector(nb_mat_int)
print(nb_mat1[,1])
print(nb_mat_int[1:9])
vars = c(0.1,0.2,0.3,0.4)
bounds_a = c(0,5)
bounds_e = c(-0.15,0.15)
bounds_t = c(0,5)
#print(as.integer(y))
dm_call_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t)

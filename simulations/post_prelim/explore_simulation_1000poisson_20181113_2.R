################################################
# Explore the simulated Poisson dataset and how robust it is
################################################

setwd('../')
source('./simulation_1000.R')
source('../analysis/post_analysis.R')
source('../multi_samples/ms_dm_call.R')

load('./simulated_neighbors_1000_poisson.RData')
load_y = TRUE
param = c(2, 0, 3)
date = '20181113'
if (!load_y) {
  y<-simulate_y_poisson_mvn(1000, sub_neighbor, param[1], param[2], param[3], 1)
  #save the current workspace with all simulated datasets
  saveRDS(y, file=paste('./post_prelim/simulated_1000poisson_',param[1],'_',param[2],'_',param[3],'_',date,'_mvn.rds', sep=""))
} else {
  y<-readRDS(file=paste('./post_prelim/simulated_1000poisson_',param[1],'_',param[2],'_', param[3],'_',date,'_mvn_2.rds', sep=""))
  #y<-readRDS(file='./simulated_1000poisson_2_0_3_ms_234.rds')
}

load_inis=FALSE
total_iter=5000
#y<-list(y1)
#when simulating a new value in each of these metropolis steps, variance of the random walk chain can be specified
#below ensures tuning of every dataset
vars<-c(0.7, 0.05, 0.09, 0.15)
n = nrow(y)
cat("data size:", n,"\n")
cat("random walk variance:", vars, "\n")
cat("Dataset II\n")
nb_mat_int = as.integer(sub_neighbor)
#there are four steps of metropolis in each iteration (including double metropolis)

#parameters for the prior distributions (uniform), for alpha, eta and tau^2
bounds_a = c(0,100)
bounds_e = c(-0.13,0.11)
bounds_t = c(0,100)
#inital guess for alpha, beta and tau^2
if (!load_inis) {
  inis1 = c(10, -0.1, 0.1)
  cat("initial values", inis1, "\n")
#saverds(inis1, file=paste('../results/',date,'/tests/','/inis_os1',counter,'.rds', sep = ''))
  winis = rnorm(n, inis1[1], sqrt(inis1[3]))
#saverds(vars, file=paste('../results/',date,'/tests/','/vars_os1_',counter,'.rds', sep = ''))
} else {
#load last results (ret object)
  tmp_winis = readRDS(paste('../results/',date,'/returned.rds', sep=''))
  ncol = length(tmp_winis$alpha)
  cat(paste("last round iteration number:", ncol, "\n"))
  winis = tmp_winis$w[,ncol]
  inis1 = c(tmp_winis$alpha[ncol], tmp_winis$eta[ncol], tmp_winis$tau2[ncol])
}

ret1 <- dm_call_wrapper(total_iter, n, y,  nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis1, winis)
jump_count = get_jump_frequency(ret1, total_iter, n)
print(jump_count)
newret <-delete_burn_in(ret1, 1000, sub_neighbor)
#cat(paste("Point estimates:\n"))
#print(mean(newret$alpha))
#print(mean(newret$eta))
#print(mean(newret$tau2))
print_param_estimates(newret)

#save ret workspace
#save(ret1, newret, vars, inis1, param, file = paste('./post_prelim/returned_workspace_',param[1], '_', param[2],'_', param[3],'_', date,'_mvn.RData', sep=""))

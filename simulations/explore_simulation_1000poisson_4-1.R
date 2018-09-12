################################################
# Explore the simulated Poisson dataset and how robust it is
################################################
#This script generates the final datasets for simulation 

setwd('./')
source('./simulation_1000.R')
source('../analysis/post_analysis.R')
source('../multi_samples/ms_dm_call.R')

load('./simulated_neighbors_1000_poisson.RData')
load_y = FALSE
param = c(2, 0,3)
n = 1000
  #net<-graph_from_adjacency_matrix(sub_neighbor, mode="undirected")
set.seed(234)
y1<-simulate_y_poisson_ms(n, sub_neighbor, param[1], param[2], param[3], 3, 5000)
load_inis=FALSE
total_iter=5000
vars<-matrix(NA, nrow=6, ncol=4)
vars[1,]<-c(0.7, 0.3, 0.15, 0.3)
cat("data size:", n,"\n")
nb_mat_int = as.integer(sub_neighbor)
#there are four steps of metropolis in each iteration (including double metropolis)
saveRDS(y1,file = './simulated_1000poisson_2_0_3_ms_234.rds')
#parameters for the prior distributions (uniform), for alpha, eta and tau^2
bounds_a = c(0,100)
bounds_e = c(-0.13,0.11)
bounds_t = c(0,100)
#inital guess for alpha, beta and tau^2
if (!load_inis) {
  inis1 = c(5,0.0,0.5)
    #saverds(inis1, file=paste('../results/',date,'/tests/','/inis_os1',counter,'.rds', sep = ''))
  winis = rnorm(n, 0, 1)
    #saverds(vars, file=paste('../results/',date,'/tests/','/vars_os1_',counter,'.rds', sep = ''))
} else {
    #load last results (ret object)
  tmp_winis = readRDS(paste('../results/',date,'/returned.rds', sep=''))
  ncol = length(tmp_winis$alpha)
  cat(paste("last round iteration number:", ncol, "\n"))
  winis = tmp_winis$w[,ncol]
  inis1 = c(tmp_winis$alpha[ncol], tmp_winis$eta[ncol], tmp_winis$tau2[ncol])
}
	
ret1 <- dm_call_wrapper(total_iter, n, y1,  nb_mat_int, vars[1,], bounds_a, bounds_e, bounds_t, inis1, winis)
jump_count = get_jump_frequency(ret1, total_iter, n)
print(jump_count)
newret <-delete_burn_in(ret1, 1000, sub_neighbor)
cat(paste("Point estimates:\n"))
print(mean(newret$alpha))
print(mean(newret$eta))
print(mean(newret$tau2))
#save ret workspace
#save(ret,vars, inis1, param, file = paste('./returned_workspace_',param[1], '_', param[2],'_', param[3],'_', code,'.RData', sep=""))

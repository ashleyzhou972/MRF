################################################
# Explore the simulated Poisson dataset and how robust it is
################################################

setwd('../')
source('./simulation_1000.R')
source('../analysis/post_analysis.R')
source('../multi_samples/ms_dm_call.R')

load('./simulated_neighbors_1000_poisson.RData')
load_y = FALSE
param = c(2, 0, 3)
date = '20181113'
if (!load_y) {
  y<-simulate_y_poisson_mvn(1000, sub_neighbor, param[1], param[2], param[3], 1)
  #save the current workspace with all simulated datasets
  saveRDS(y, file=paste('./post_prelim/simulated_1000poisson_',param[1],'_',param[2],'_',param[3],'_',date,'_mvn_2.rds', sep=""))
} else {
  #y<-readRDS(file=paste('./post_prelim/simulated_1000poisson_',param[1],'_',param[2],'_', param[3],'_',date,'_mvn.rds', sep=""))
  y<-readRDS(file='./simulated_1000poisson_2_0_3_ms_234.rds')
}

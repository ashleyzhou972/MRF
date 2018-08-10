###################
# Comparison between models 
# 500 hic threshold vs 5000 hic threshold
# results for hic 500 is in 20180727  (1000 iterations)
# results for hic 5000 is in 20180731 (2000 iterations)
##############################

#setwd('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/analysis')
library(parallel)
library(Matrix)
source('./post_analysis.R')
####read results 1
ret500<-readRDS('ret500.rds')

neighbor500<-readRDS('neighbor500.rds')

ret5000<-readRDS('ret5000.rds')

neighbor5000<-readRDS('neighbor5000.rds')


new_ret500<-delete_burn_in(ret500,170,as.matrix(neighbor500))
new_ret5000<-delete_burn_in(ret5000,200, as.matrix(neighbor5000))

N=19452

newret = list(new_ret500, new_ret5000)
no_cores = 15
cl <-makeCluster(no_cores)
yhat<-parLapply(cl, newret, fun=get_fitted_y, N=N, iters=5)
#yhat<-lapply(newret, FUN=get_fitted_y, N=N, iters=2)

yhat500 = yhat[[1]]
yhat5000 = yhat[[2]]

saveRDS(yhat500, file="./yhat500.rds")
saveRDS(yhat5000, file="./yhat5000.rds")




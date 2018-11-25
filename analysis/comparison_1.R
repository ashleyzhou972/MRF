###################
# Comparison between models 
# 500 hic threshold vs 5000 hic threshold
# results for hic 500 is in 20180727  (1000 iterations)
# results for hic 5000 is in 20180731 (2000 iterations)
##############################

#setwd('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/analysis')

source('./post_analysis.R')
####read results 1
date = '20180727'
#load(paste('../results/',date,'/returned.RData', sep=""))
#save as rds file for ease in loading
#saveRDS(ret1, file='ret500.rds')
ret500<-readRDS('ret500.rds')

thres = 500
#load neighbor
#load(paste('/home/nzhou/hic/rao2014/adj_mat_',thres,'.RData', sep=""), verbose=T)
#saveRDS(adj_mat_500, file="neighbor500.rds")
neighbor500<-readRDS('neighbor500.rds')

####read results 2
date = '20180731'
#load(paste('../results/',date,'/returned.RData', sep=""))
#save as rds file 
#saveRDS(ret1, file='ret5000.rds')
ret5000<-readRDS('ret5000.rds')

thres=5000
#load neighbor
#load(paste('/home/nzhou/hic/rao2014/adj_mat_',thres,'.RData', sep=""), verbose=T)
#saveRDS(adj_mat_5000, file="neighbor5000.rds")
neighbor5000<-readRDS('neighbor5000.rds')


new_ret500<-delete_burn_in(ret500,170)
new_ret5000<-delete_burn_in(ret5000,200)

N=19452

yhat500<-get_fitted_y(N, mean(new_ret500$alpha), mean(new_ret500$eta), neighbor500)
yhat5000<-get_fitted_y(N, mean(new_ret5000$alpha), mean(new_ret5000$eta), neighbor5000)

diff = yhat500-yhat5000
order(diff, decreasing=T)



#comparing two neighborhoods

net500<- graph_from_adjacency_matrix(neighbor500)
net5000<- graph_from_adjacency_matrix(neighbor5000)

ecount(net500)
ecount(net5000)


par(mfrow=c(3,1))
hist(new_ret500$alpha,prob=T, xlab = "alpha", main="")
hist(new_ret500$eta,prob=T, xlab="eta", main="")
hist(new_ret500$tau2,prob=T,xlab="tau2", main="")


hist(new_ret5000$alpha,prob=T, xlab = "alpha", main="")
hist(new_ret5000$eta,prob=T, xlab="eta", main="")
hist(new_ret5000$tau2,prob=T,xlab="tau2", main="")

plot(ret500$eta, type='l', xlab="iteration", ylab="eta")
abline(v=170, col=2)
abline(h = mean(new_ret500$eta), col=3)
plot(ret5000$eta, type='l', xlab="iteration", ylab="eta")
abline(v=200, col=2)
abline(h = mean(new_ret5000$eta), col=3)

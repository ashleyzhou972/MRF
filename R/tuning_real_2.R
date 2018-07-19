########
#using Rao 2014 data
#threshold = 500
########


ptm<-proc.time()

setwd('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/R')
source('../../read.R')
source('../c/dm_call.R')
source('../../functions_for_metropolis.R')
source('../../scaleReductionFactor.R')
source('./simulation_1000.R')
#source('./pseudolikelihood.R')

date = '20180718'
thres = 500
load(paste('../../../../rao2014/adj_mat_',thres,'.RData', sep = ''))
c<-adj_mat_500
cat("Data loaded\n")

#cat("Calculating eigen values\n")
#evalues<-eigen(c)$value
#save(evalues, file = '../results/20180626/evalues.RData')
#eta_max <- 1/max(evalues)
#eta_min <- 1/min(evalues)
eta_max = 0.010709986
eta_min = -0.013782385
#cat("Eigen value calculation complete\n")
cat(paste('max eta:' ,eta_max),'\n')
cat(paste('min eta:' , eta_min),'\n')

sub_neighbor<-as.matrix(c)
y<-data[data$Ensembl_ID%in%colnames(c),"count"]


#see degree distribution
sub_net = graph_from_adjacency_matrix(sub_neighbor,mode='undirected')
#barplot(table(degree(sub_net)), xlab='Degree', ylab="Frequency",  col="#C8102E", main = pathway_name)
#hist(degree(sub_net), breaks = 100, xlab = "Degree", ylim = c(0,4000), col="#C8102E",main='all',cex.axis = 1.5 , cex.lab=1.5, cex.main = 2)


total_iter = 50
B = 1000


N = length(y)
save(N,file=paste('../results/',date,'/N.RData', sep = ''))
nb_mat_int = as.integer(sub_neighbor)


#limits for eta: -0.079, 0.022
#There are four steps of metropolis in each iteration (including double metropolis)
#When simulating a new value in each of these metropolis steps, variance of the random walk chain can be specified
vars = c(0.8,0.07,0.007,0.1)
save(vars, file=paste('../results/',date,'/vars.RData', sep = ''))
#first var is ignored in the gaussian case
#parameters for the prior distributions (uniform), for alpha, eta and tau^2
bounds_a = c(0,10)
bounds_e = c(eta_min,eta_max)
bounds_t = c(0,10)
#inital guess for alpha, beta and tau^2

inis1 = c(5,0.0,0.5)
ret1 <- dm_call_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis1)
jump_count = get_jump_frequency(ret1, total_iter, N)
save(ret1, file = paste('../results/',date,'/returned.RData', sep = ''))
save(jump_count, file = paste('../results/',date,'/jumps.RData', sep = ''))

proc.time()-ptm
duration=proc.time()-ptm
cat('Time elapsed',duration[3],'\n')
save(duration, file = paste('../results/',date,'/duration.RData', sep = ''))



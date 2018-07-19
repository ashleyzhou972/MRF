########
#using Jin 2013 data
########
ptm<-proc.time()

setwd('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/R')
source('../../read.R')
source('../c/dm_call.R')
source('../../functions_for_metropolis.R')
source('../../scaleReductionFactor.R')
source('./simulation_1000.R')
#source('./pseudolikelihood.R')


y<-data$count
sub_neighbor<-c

cat("Data loaded\n")
total_iter = 3000
B = 1000


N = length(y)
nb_mat_int = as.integer(sub_neighbor)

#limits for eta: -0.079, 0.022
#There are four steps of metropolis in each iteration (including double metropolis)
#When simulating a new value in each of these metropolis steps, variance of the random walk chain can be specified
vars = c(0.8,0.1,7.0,0.1)
save(vars, file='../results/20180622/vars.RData')
#first var is ignored in the gaussian case
#parameters for the prior distributions (uniform), for alpha, eta and tau^2
bounds_a = c(0,10)
bounds_e = c(-0.08,0.022)
bounds_t = c(0,10)
#inital guess for alpha, beta and tau^2

inis1 = c(5,0.0,0.5)
ret1 <- dm_call_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis1)
jump_count = get_jump_frequency(ret1, total_iter, N)
save(ret1, file = '../results/20180622/returned.RData')
save(jump_count, file = '../results/20180622/jumps.RData')

proc.time()-ptm
duration=proc.time()-ptm
cat('Time elapsed',duration[3],'\n')
save(duration, file = '../results/20180622/duration.RData')



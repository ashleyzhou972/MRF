setwd('/home/nzhou/hic/IMR90/work/dm/MRF_HIC_GE/c')
source('./dm_call.R')
source('../../../functions_for_metropolis.R')
#Total number of iterations
total_iter = 10000
B = 2000
#load observed data
load('./simulated_y__new.RData')
load('./simulated_neighbors_new.RData')
#True parameter values
# alpha = 2
# eta = 0.1
# tau2 = 2
N = length(y)
nb_mat_int = as.integer(sub_neighbor)
#There are four steps of metropolis in each iteration (including double metropolis)
#When simulating a new value in each of these metropolis steps, variance of the random walk chain can be specified
vars = c(0.6,0.5,0.03,0.14)
#parameters for the prior distributions (uniform), for alpha, eta and tau^2
bounds_a = c(0,10)
bounds_e = c(-0.20,0.155)
bounds_t = c(0,10)
#inital guess for alpha, beta and tau^2
inis = c(2,0.1,1)

ptm<-proc.time()

ret <- dm_call_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis)
#results
alpha = ret[[2]]
eta = ret[[3]]
tau2 = ret[[4]]
jump_frequency_w = ret[[5]][1]/(N*total_iter)
jump_frequency_alpha = ret[[5]][2]/total_iter
jump_frequency_eta = ret[[5]][3]/total_iter
jump_frequency_tau2 = ret[[5]][4]/total_iter
jump_count = c(jump_frequency_w, jump_frequency_alpha, jump_frequency_eta, jump_frequency_tau2)
cat('Jump rate for w is ', jump_frequency_w, '\n')
cat('Jump rate for alpha is ', jump_frequency_alpha, '\n')
cat('Jump rate for eta is ', jump_frequency_eta, '\n')
cat('Jump rate for tau2 is ', jump_frequency_tau2, '\n')

proc.time()-ptm
duration=proc.time()-ptm
cat('Time elapsed',duration[3],'\n')
#dm for regular metropolis
subDir = paste('../../../mcmc_results/dm/',format(Sys.Date(),"%Y%m%d"),sep='')
if (!dir.exists(subDir))  dir.create(subDir)
setwd(subDir)
save(alpha,file = './sim_alpha.RData') 
save(eta,file = './sim_eta.RData')
save(tau2,file = './sim_tau2.RData')
save(jump_count,file='./sim_jumpcount.RData')

png(filename = './sim_alpha_iterations_plot.png')
plot_iterations(alpha,alpha0,0,20,'alpha')
dev.off()
png(filename = './sim_eta_iterations_plot.png')
plot_iterations(eta,eta0,-0.2,0.2,'eta')
dev.off()
png(filename = './sim_tau2_iterations_plot.png')
plot_iterations(tau2,tau20,0,21,'tau2')
dev.off()

alpah_post = get_posterior_distribution(B,alpha)
eta_post = get_posterior_distribution(B,eta)
tau2_post = get_posterior_distribution(B,tau2)

png(filename='./sim_eta_hist.png')
hist(eta_post)
dev.off()

print(summary(eta_post))
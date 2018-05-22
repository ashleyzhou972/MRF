setwd('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/R')
source('../c/dm_call.R')
source('../../functions_for_metropolis.R')
source('../../scaleReductionFactor.R')
#Total number of iterations
total_iter = 10000
B = 2000
#load observed data
load('./simulated_y__new.RData')
load('./simulated_neighbors_new.RData')
#True parameter values
alpha_true = 2
eta_true = 0.1
tau2_true = 2
N = length(y)
nb_mat_int = as.integer(sub_neighbor)
#There are four steps of metropolis in each iteration (including double metropolis)
#When simulating a new value in each of these metropolis steps, variance of the random walk chain can be specified
vars = c(0.65,0.3,0.02,0.13)
#parameters for the prior distributions (uniform), for alpha, eta and tau^2
bounds_a = c(0,10)
bounds_e = c(-0.20,0.155)
bounds_t = c(0,10)
#inital guess for alpha, beta and tau^2

inis1 = c(0.5,-0.16,0.5)
inis2 = c(2,0.05,2)
inis3 = c(5,0,5)
inis4 = c(10,-0.1,10)
inis5 = c(1,0.19,1)

ptm<-proc.time()


ret1 <- dm_call_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis1)
ret2 <- dm_call_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis2)
ret3 <- dm_call_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis3)
ret4 <- dm_call_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis4)
ret5 <- dm_call_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis5)


#results
alpha = cbind(ret1[[2]],ret2[[2]],ret3[[2]],ret4[[2]],ret5[[2]])
eta = cbind(ret1[[3]],ret2[[3]],ret3[[3]],ret4[[3]],ret5[[3]])
tau2 = cbind(ret1[[4]],ret2[[4]],ret3[[4]],ret4[[4]],ret5[[4]])

#scale reduction factor
scaleReductionFactor(alpha)
scaleReductionFactor(eta)
scaleReductionFactor(tau2)

#jump rate
jump_count = get_jump_frequency(ret1, total_iter, N)

#plot
plot_iterations(total_iter, ret1$alpha, inis1[1],bounds_a[1], bounds_a[2],"alpha")
plot_iterations(total_iter, ret1$eta, inis1[2],bounds_e[1], bounds_e[2],"eta")
plot_iterations(total_iter, ret1$tau2, inis1[3],bounds_t[1], bounds_t[2],"tau2")


par(mfrow = c(3,1))
alpha_main = paste("alpha (true= ", alpha_true, ")", sep = "")
plot_iterations(total_iter, ret1$alpha, inis1[1],bounds_a[1], bounds_a[2],"alpha", alpha_main, alpha_true)
eta_main = paste("eta (true= ", eta_true, ")", sep = "")
plot_iterations(total_iter, ret1$eta, inis1[2],bounds_e[1], bounds_e[2],"eta", eta_main, eta_true)
tau2_main = paste("tau2 (true= ", tau2_true, ")", sep = "")
plot_iterations(total_iter, ret1$tau2, inis1[3],bounds_t[1], bounds_t[2],"tau2", tau2_main, tau2_true)

################################################################################
#Below for running Rscript in command line
################################################################################


for (i in 1:5) {
  cat('Jump rate for w for chain ', i, ' is ', jump_frequency_w[i], '\n')
  cat('Jump rate for alpha for chain ', i, ' is ', jump_frequency_alpha[i], '\n')
  cat('Jump rate for eta for chain ', i, ' is ', jump_frequency_eta[i], '\n')
  cat('Jump rate for tau2 for chain ', i,  ' is ', jump_frequency_tau2[i], '\n')
}

proc.time()-ptm
duration=proc.time()-ptm
cat('Time elapsed',duration[3],'\n')
#dm for regular metropolis
subDir = paste('../../mcmc_results/dm/',format(Sys.Date(),"%Y%m%d"),sep='')
if (!dir.exists(subDir))  dir.create(subDir)
setwd(subDir)
alpha = ret1$alpha
save(alpha,file = './sim_alpha.RData') 
eta = ret1$eta
save(eta,file = './sim_eta.RData')
tau2 = ret1$tau2
save(tau2,file = './sim_tau2.RData')
save(jump_count,file='./sim_jumpcount.RData')

alpha0 = c(inis1[1], inis2[1], inis3[1], inis4[1], inis5[1])
eta0 = c(inis1[2], inis2[2], inis3[2], inis4[2], inis5[2])
tau20 = c(inis1[3], inis2[3], inis3[3], inis4[3], inis5[3])
png(filename = './sim_alpha_iterations_plot.png')
plot_iterations_srf(total_iter, alpha,alpha0,0,bounds_a[2],'alpha')
dev.off()
png(filename = './sim_eta_iterations_plot.png')
plot_iterations_srf(total_iter, eta,eta0,bounds_e[1],bounds_e[2],'eta')
dev.off()
png(filename = './sim_tau2_iterations_plot.png')
plot_iterations_srf(total_iter, tau2,tau20,0,bounds_t[2],'tau2')
dev.off()

alpah_post = get_posterior_distribution_srf(total_iter, B,alpha)
eta_post = get_posterior_distribution(total_iter, B,eta)
tau2_post = get_posterior_distribution(total_iter, B,tau2)

png(filename='./sim_eta_hist.png')
hist(eta_post)
dev.off()

print(summary(eta_post))

##########################################################
# Simulated Datasets run on multiple chains
########################################################

setwd('./')
source('../multi_samples/ms_dm_call.R')
source('../R/functions_for_metropolis.R')
source('./simulation_1000.R')
#source('./pseudolikelihood.R')
#Total number of iterations
total_iter = 5000
true_params <-matrix(c(1,-0.05,1,1,0.05,1,2,-0.1,3,2,0.1,3,2,0,3), byrow = T, ncol=3)
#Round is the round of MCMC
#Since MCMC could be continued using last round's values as inital values
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Must input two arguments")
} else {
  round = as.numeric(args[1])
  true_param = true_params[as.numeric(args[2]),]
}
load('./simulated_neighbors_1000_poisson.RData')
y = readRDS(paste('./simulated_1000poisson_', true_param[1],'_', true_param[2], '_', true_param[3],'_ms_234.rds', sep=''))
vars = matrix(NA, nrow = 5, ncol = 4)
vars[1,]<-c(0.7, 0.1, 0.15, 0.1)
vars[2,]<-c(0.7, 0.3, 0.15, 0.25)
vars[3,]<-c(0.7, 0.3, 0.15, 0.25)
vars[4,]<-c(0.7, 0.15, 0.17, 0.05)
vars[5,]<-c(0.7, 0.15, 0.17, 0.05) 

N = nrow(y)
nb_mat_int = as.integer(sub_neighbor)
bounds_a = c(0,100)
bounds_e = c(-0.13,0.11)
bounds_t = c(0,100)
#inital guess for alpha, beta and tau^2
if (round==1) {
  load_inis = FALSE
} else {
  load_inis = TRUE
}

if (!load_inis) {
  inis1 = c(5,0.0,5)
  inis2 = c(0,0.0,0)
  inis3 = c(10,0.2,10)
  inis4 = c(2, -0.2, 2)
  inis5 = c(2, 0.2, 2)
  inis_all = rbind(inis1,inis2, inis3, inis4, inis5)
  winis_one = rnorm(N, 0, 1)
  #repeat the same w initial for all five chains
  winis_all = matrix(winis_one, nrow = 5, ncol = N, byrow = TRUE)
}  else {
  #load last results (ret object)
  load(paste("./results_multiple_chains_round",round-1, "_" ,true_param[1], "_", true_param[2],"_", true_param[3],".RData", sep=""))
  old_ret <- ret
  inis_all = matrix(NA, nrow = 5, ncol =3)
  winis_all = matrix(NA, nrow = 5, ncol = N)
  for (j in 1:5) {
    alpha_pre = ret[[j]]$alpha
    eta_pre = ret[[j]]$eta
    tau2_pre = ret[[j]]$tau2
    ncol = length(alpha_pre)
    winis = ret[[j]]$w[,ncol]
    inis = c(alpha_pre[ncol], eta_pre[ncol], tau2_pre[ncol])
    cat("Initial parameters are ", inis, "\n")
    inis_all[j,]<-inis
    winis_all[j,]<-winis
  }
  remove("ret")
}

cat("Arguments finished loading.\n")
new_ret <-list()
new_jumps<-list()
for (i in 1:5) {
  cat(paste("Now running chain", i, "\n"))
  ret1 <- dm_call_wrapper(total_iter, N, y,  nb_mat_int, vars[i,], bounds_a, bounds_e, bounds_t, inis_all[i,], winis_all[i,])
  jump_count = get_jump_frequency(ret1, total_iter, N)
  new_ret[[i]] <- ret1
  new_jumps[[i]]<-jump_count
  print(jump_count)
}

if (round==1) {
  ret = new_ret
}  else {
  ret <-list()
  for (j in 1:5) {
    ret[[j]] = assemble_rets(old_ret[[j]], new_ret[[j]])
  }
}

save(ret, file = paste("./results_multiple_chains_round",round, "_" ,true_param[1], "_", true_param[2],"_", true_param[3],".RData", sep=""))

#calculate scale reduction factor
SRF_alpha_input = cbind(ret[[1]]$alpha, ret[[2]]$alpha, ret[[3]]$alpha, ret[[4]]$alpha, ret[[5]]$alpha)
SRF_eta_input = cbind(ret[[1]]$eta, ret[[2]]$eta, ret[[3]]$eta, ret[[4]]$eta, ret[[5]]$eta)
SRF_tau2_input = cbind(ret[[1]]$tau2, ret[[2]]$tau2, ret[[3]]$tau2, ret[[4]]$tau2, ret[[5]]$tau2)
srf_alpha = scaleReductionFactor(SRF_alpha_input)
srf_eta = scaleReductionFactor(SRF_eta_input)
srf_tau2 = scaleReductionFactor(SRF_tau2_input)
cat("Scale Reduction Factors:\n")
cat(srf_alpha,"\n")
cat(srf_eta,"\n")
cat(srf_tau2,"\n")

if (srf_alpha<1.1 & srf_eta<1.1 & srf_tau2 <1.1) {
  cat("Convergence reached\n")
  cat(TRUE, "\n")
}  else {
  cat(FALSE, "\n")
}

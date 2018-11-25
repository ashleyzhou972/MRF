###############################
#Run 100 random subsets (node subsets) and compare \hat{eta}
#number of nodes in the random subsets is not fixed
###############################
setwd('./')
source('../R/read.R')
source('../multi_samples/ms_dm_call.R')
source('./post_analysis.R')

library(igraph)
neighbor5000<-readRDS('neighbor5000.rds')
net5000<- graph_from_adjacency_matrix(neighbor5000, mode="undirected")
#signal pathway
glist<-as.character(read.table('../../rao2014/pathway/genes_signal_GO.txt', sep="\n", header=F)$V1)
#Because this HiC data doesn't contain the Y chromosome
#Need to filter this gene list
signal = colnames(neighbor5000)[colnames(neighbor5000)%in%glist]

total_iter = 5000
B = 1

bounds_a<-c(0,100)
bounds_t<-c(0,100)
inis_par<-c(1,0,1)

#######################
y_ob<-data[data$Ensembl_ID%in%signal,'count']
neighbor_ob<-as.matrix(neighbor5000[which(colnames(neighbor5000)%in%signal), which(colnames(neighbor5000)%in%signal)])
evalues_ob<-eigen(neighbor_ob)$value
eta_min_ob = 1/min(evalues_ob)
eta_max_ob = 1/max(evalues_ob)
cat("eta bounds for observed are", eta_min_ob, eta_max_ob, "\n")
bounds_e_ob = c(eta_min_ob, eta_max_ob)
inis_w_ob<-rnorm(length(signal), 0, 1)

vars_ob = c(0.7,0.2,(eta_max_ob-eta_min_ob)/2-0.08,0.2)

ret_ob<-dm_call_wrapper(total_iter, length(signal), y_ob, as.integer(neighbor_ob), vars_ob, bounds_a, bounds_e_ob, bounds_t, inis_par, inis_w_ob)
jumps_ob<-get_jump_frequency(ret_ob, total_iter, length(signal))
print(jumps_ob)
newret_ob<-delete_burn_in(ret_ob, B, neighbor_ob)
eta_l_ob = t.test(newret_ob$eta)$conf.int[1]
eta_u_ob = t.test(newret_ob$eta)$conf.int[2]
cat("Estimated eta from observed is", mean(newret_ob$eta), "(", eta_l_ob, ",", eta_u_ob, ")\n")
save(ret_ob, eta_min_ob, eta_max_ob, jumps_ob, file="./ret_observed_signal.RData")
saveRDS(ret_ob, file = "./ret_observed_signal.RDS")
##########################33


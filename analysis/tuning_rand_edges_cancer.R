###############################
#Run 100 random subsets (edge subsets) and compare \hat{eta}
#number of nodes in the random subset is set
###############################
setwd('./')
source('./random_subsets.R')
source('../R/read.R')
source('../multi_samples/ms_dm_call.R')
source('./post_analysis.R')

#args = commandArgs(trailingOnly=T)
#mode = args[1]
mode = "cancer"

neighbor5000<-readRDS('neighbor5000.rds')
net5000<- graph_from_adjacency_matrix(neighbor5000, mode="undirected")
if (mode=="cancer") {
  glist<-as.character(read.table('../../rao2014/pathway/genes_cancer_pathway', sep="\n", header=F)$V1) 
}  else if (mode =="signal") {
  glist<-as.character(read.table('../../rao2014/pathway/genes_signal_GO.txt', sep="\n", header=F)$V1)
}
cancer<-colnames(neighbor5000)[colnames(neighbor5000)%in%glist]

total_iter = 3000
B = 200
N = 1000 #100 random subsets
rand_subsets=make_subset_edges_gnm(cancer, neighbor5000,N)

bounds_a<-c(0,100)
bounds_t<-c(0,100)
inis_par<-c(1,0,1)
n = length(cancer)
#inis_w<-rnorm(n, 0, 1)
#saveRDS(inis_w, file = paste("./inital_w_",n,"_edges.rds", sep=""))
inis_w<-readRDS(paste("./inital_w_",n,"_edges.rds", sep=""))

all_eta<-rep(NA, N)
all_eta_CI<-matrix(NA, nrow=N, ncol=5)
#all_ret<-list()
i = 1
while(i <= N) {
  nb<-rand_subsets[[i]]
  neighbor_int<-as.integer(nb)
  nodes<-colnames(nb)
  y<-data[data$Ensembl_ID%in%nodes,'count']
  
  evalues<-eigen(nb)$value
  eta_min = 1/min(evalues)
  eta_max = 1/max(evalues)
  
  vars = c(0.8,0.6,(eta_max-eta_min)/2,0.5)
  bounds_e<-c(eta_min, eta_max)
  
  ret<-dm_call_wrapper(total_iter,n, y, neighbor_int, vars, bounds_a, bounds_e, bounds_t, inis_par, inis_w)
  jumps<-get_jump_frequency(ret, total_iter, n)
  jumps_vec<-c(jumps$w, jumps$alpha, jumps$eta, jumps$tau2)
  #all_ret[[i]]<-ret
  if (sum(jumps_vec<0.7)==4 & sum(jumps_vec>0.2)==4) {
    newret<-delete_burn_in(ret, B, nb)
    eta_hat<-mean(newret$eta)
    eta_ttest<-t.test(newret$eta)
    eta_l<-round(eta_ttest$conf.int[1],5)
    eta_u<-round(eta_ttest$conf.int[2],5)
    cat("Estimated eta is", mean(newret$eta), "(", eta_l, ",", eta_u, ")\n")
    all_eta[i]<-eta_hat
    #all_ret[[i]]<-ret
    all_eta_CI[i,]<-c(eta_hat, eta_l, eta_u, eta_min, eta_max)
    i = i+1
  } else {
    cat("Jump counts are invalid, results are disregarded\n")
  } 
  rm("ret")
}
###################################################3
save.image(file=paste("./rand_edges_",mode,"_workspace.RData", sep=""))
saveRDS(all_eta_CI, file = paste("./rand_edges_",mode,"_all_eta_CI", sep=""))


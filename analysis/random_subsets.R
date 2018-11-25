#############################################
# Random subsets of data (sub neighborhoods)
# Two types of random subsets of data
  # - For a list of pathway genes, generate random edges (random graph)
  # - For the observed network structure, sample (stratified) a random list of nodes (same number as pathway genes)
#############################################

library(igraph)

#source('./post_analysis.R')
#neighbor5000<-readRDS('neighbor5000.rds')
#net5000<- graph_from_adjacency_matrix(neighbor5000, mode="undirected")
#cancer<-as.character(read.table('/home/nzhou/hic/rao2014/pathway/genes_cancer_pathway', sep="\n", header=F)$V1)


#random sampling of the nodes

#stratified sampling of the nodes
#strata by degree:
#0-19, 20-39, 40-59, 60-79, 80-Inf
make_strata<-function(net){
  degrees=degree(net)
  br = seq(min(degrees), max(degrees)+20, 20)
  histo=hist(degrees, breaks=br, plot=F)
  strata<-cbind(br[1:(length(br)-1)], br[-1])
  colnames(strata)<-c("stratum_start", "stratum_end")
  return(strata)
}

#################################################################3

sample_subset_nodes_statified<-function(strata, num_strata,degrees, prop){
  subset = c() #subset is just 
  for (stratum in 1:num_strata){
    node_ids = which(degrees>=strata[stratum,1] & degrees<strata[stratum, 2])
    count = length(node_ids)
    sampled = sample(node_ids,floor(prop*count))
    subset = c(subset, sampled)
  }
  
  return(subset)
}

make_subset_nodes_statified<-function(neighbor, strata,n, N) {
  #Do the sampling N times
  prop = n/nrow(neighbor)
  net<-graph_from_adjacency_matrix(neighbor, mode="undirected")
  degrees = degree(net)
  num_strata = nrow(strata)
  all_subnets<-list()
  for (i in 1:N) {
    sublist <- sample_subset_nodes_statified(strata, num_strata, degrees, prop)  
    subnet<-as.matrix(neighbor[which(colnames(neighbor)%in%names(sublist)),which(colnames(neighbor)%in%names(sublist))])
    all_subnets[[i]]<-subnet
  }
  return(all_subnets)
}

############################################################
sample_subset_edges_gnm<-function(n,m) {
  tmp_net<-sample_gnm(n,m, directed=F)
  return(tmp_net)
}
make_subset_edges_gnm<-function(glist, neighbor, N) {
  n = length(glist)
  sub_neighbor<-as.matrix(neighbor[which(colnames(neighbor)%in%glist),which(colnames(neighbor)%in%glist)])
  subnet = graph_from_adjacency_matrix(sub_neighbor, mode = "undirected")
  m = ecount(subnet)
  all_subnets<-list()
  for (i in 1:N) {
    tmp_net <-sample_subset_edges_gnm(n,m)
    sub_net<-set_vertex_attr(tmp_net, "name", value = as.character(glist))
    all_subnets[[i]]<-as.matrix(as_adjacency_matrix(sub_net, type="both", names = T))
  }
  return(all_subnets)
}



######################
# Test
######################

#Type2
#strata<-make_strata(net5000)
#n = length(cancer)
#t1<-make_subset_nodes_statified(neighbor5000, strata = strata, n , 2)


#Type1
#t2<-make_subset_edges_gnm(cancer, neighbor5000, 2)









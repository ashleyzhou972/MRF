###################################################
#This script performs the same MRF analysis 
#instead of with the full 19676 genes
#only select genes in a single pathway
#updated 20180622
###################################################
#library(igraph, warn.conflicts = F, quietly = T)


ptm<-proc.time()

setwd('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/R')
#read.R no longer loads neighborhood matrix
source('../../read.R')
source('../c/dm_call.R')
#source('../../functions_for_metropolis.R')
source('../../scaleReductionFactor.R')
source('./simulation_1000.R')
source('./pseudolikelihood.R')

#load neighborhood matrix by threshold
thres = 500
load(paste('../../../../rao2014/adj_mat_',thres,'.RData', sep = ''))
c<-adj_mat_500

#load neighborhood matrix from Jin et al
#load('../../c_matrix.RData')

#read in gene list in pathway

date = '20180626'
#pathway_names = c('cancer', 'smooth','TGF','TNF')
#pathway_names = c('smooth')
#pathway_names = c('cancer')
pathway_names = c('TGF')
#pathway_names = c('TNF')

total_iter = 10000
#B = 1000


for (pathway_name in pathway_names){
	cat(paste('Now runnig pathway', pathway_name,'\n'))
	glist=as.vector(read.table(paste('../../../../rao2014/pathway/genes_',pathway_name,'_pathway', sep='')))
	#map and filter
	#check colnames and rownames are the same
	#sum(colnames(c) == rownames(c))
	geneList = colnames(c)[colnames(c)%in%glist[,1]]
	y<-data[data$Ensembl_ID%in%geneList,'count']
	sub_neighbor<-as.matrix(c[which(colnames(c)%in%glist[,1]),which(colnames(c)%in%glist[,1])])

	#read adjacency matrix as graph and check edge
	sub_net = graph_from_adjacency_matrix(sub_neighbor, mode = 'undirected')
	png(paste('/home/nzhou/hic/presentations/degree_',pathway_name,'.png', sep=''))
	barplot(table(degree(sub_net)), xlab='Degree', ylab="Frequency",  col="#C8102E", main = pathway_name, 
  	        cex.axis = 1.5 , cex.names=1.5, cex.lab=1.5, cex.main = 2, ylim = c(0,80))
	dev.off()
}
	cat(paste('vcount:', vcount(sub_net),'\n'))
	cat(paste('ecount:', ecount(sub_net),'\n'))
	evalues<-eigen(sub_neighbor)$value
	eta_max <- 1/max(evalues)
	eta_min <- 1/min(evalues)
	cat(paste('max eta:' ,eta_max),'\n')
	cat(paste('min eta:' , eta_min),'\n')

	cat("Data loaded\n")
	
	N = length(y)
	nb_mat_int = as.integer(sub_neighbor)

	#There are four steps of metropolis in each iteration (including double metropolis)
	#When simulating a new value in each of these metropolis steps, variance of the random walk chain can be specified
	vars = c(0.8,0.6,0.2,0.5)
	foldername = paste('../results/',pathway_name,'/',date,'/', sep='')
	save(vars, file=paste(foldername,'vars.RData', sep = ''))
	#first var is ignored in the gaussian case
	#parameters for the prior distributions (uniform), for alpha, eta and tau^2
	bounds_a = c(0,10)
	bounds_e = c(eta_min,eta_max)
	bounds_t = c(0,10)
	#inital guess for alpha, beta and tau^2

	inis1 = c(0,0.0,0.5)
	save(inis1, file=paste(foldername, 'inis.RData',sep=''))
	ret1 <- dm_call_wrapper(total_iter, y, nb_mat_int, vars, bounds_a, bounds_e, bounds_t, inis1)
	jump_count = get_jump_frequency(ret1, total_iter, N)
	save(ret1, file = paste(foldername,'returned.RData', sep = ''))
	save(jump_count, file = paste(foldername,'jumps.RData', sep = ''))

	proc.time()-ptm
	duration=proc.time()-ptm
	cat('Time elapsed',duration[3],'\n')
	save(duration, file = paste(foldername,'duration.RData', sep = ''))
}


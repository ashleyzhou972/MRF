#read data
read_file <- function(filename){
  tsv = read.table(filename,header = T, sep = "\t")
  return(tsv)
}

#clean data
clean_data<-function(count_vector){
  #round up non-integer values
  ids = which(count_vector%%1!=0)
  count_vector[ids] = round(count_vector[ids])
  return(count_vector)
}
raw = read_file('./rnaseq_corrected.txt')
count = clean_data(raw[,2])
data = cbind(raw[,-2],count)

## Only 16196 genes are mapped to an interacting loci
# genes = read.table(sep='\n',file = "./genes_mapped.txt")
# ensembl = genes[,1]
# newdata = data[data$Ensembl_ID%in% ensembl,]
# data = newdata

#inspect
nonzero_t = data[data$count>=1,]
log_count = log(nonzero_t$count)
nonzero = cbind(nonzero_t,log_count)
rm(nonzero_t)
zerocount = sum(data$count==0)
#hist(nonzero$log_count)
#hist(data[data$housekeeping==1,'count'])

#neighborhood matrix
#load('/home/nzhou/hic/IMR90/work/c_matrix.RData')

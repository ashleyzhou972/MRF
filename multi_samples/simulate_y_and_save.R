library(igraph)
source('../R/simulation_1000.R')
truevalues = matrix(c(2,-0.1,4,3,0.1,2,1,0.115,7,1,-0.13,2,2,0,5), nrow = 5, byrow = T)
load('../R/simulated_neighbors_1000_poisson.RData')
net = graph_from_adjacency_matrix(sub_neighbor, mode = "undirected")
set.seed(4)
for (counter in 5:5){
  alpha_true= truevalues[counter,1]
  eta_true = truevalues[counter,2]
  tau2_true = truevalues[counter,3]
  filename = paste("../R/new_simulated_1000poisson_ms_",alpha_true,"_", eta_true, "_", tau2_true,".rds",sep = "")
  y=simulate_y_poisson_ms(net,alpha_true,eta_true,tau2_true,3,10000)
  saveRDS(y,file=filename)
  cat("Simulation done and saved\n")
}

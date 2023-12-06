# -------------------------------------------------------------------
# Script: Alternative Bootstrap Methods
# -------------------------------------------------------------------
# Description:
# This script applies a range of alternative bootstrapping methods to generated 
# networks. It is designed to evaluate and compare the effectiveness of these
# methods in estimating the standard errors of certain graph statistics.
#
# Methods Implemented:
# A. Empirical Graphon Bootstrap - Applies bootstrap techniques based on empirical graphon estimates.
# B. Subsampling Method - Utilizes subsampling for network analysis to handle large-scale networks.
# C. Embedding Method - Employs embedding techniques for dimensional reduction and network visualization.
# D. Neighborhood-Sampling Method - Focuses on sampling methods within neighborhood structures of networks.
# E. Proposed Method with Edge Sampling - A novel method emphasizing edge sampling in network analysis.
# F. Proposed Method with Random Neighbors - Introduces a new approach using random neighbor selection.
# G. Histogram Bootstrap - Implements a histogram-based bootstrap approach for network data.
#
#
# Author: Tianhai Zu
# Affiliation: University of Texas at San Antonio
# Created on: 12/05/2023
#
# Requirements:
# The script requires the following R packages: 
# - igraph (for network generation and analysis)
# - parallel, foreach (for parallel computation)
# - localboot (this package)
# Install these packages using install.packages() if not already installed.
# Additional files, like real-world network data, should be in the
# specified directory.

# Library Imports
library(igraph)
library(parallel)
library(foreach)
#library(localboot) #currently commented out

# temp source, will be deleted
source('../graph_utils.r')
source('../graph_boot_funs.R')

# -------------------------------------------------------------------
# Script Configuration
# -------------------------------------------------------------------
# Define any global variables, parameters, or method-specific settings
# here. This section should include configurations for each of the
# alternative methods being applied.

# -------------------------------------------------------------------
# Script Starts Here
# -------------------------------------------------------------------

#hist 
size = 200
M = 48
#M_true = 2400
B = 500
node_sampling = TRUE
cl_nodes = 10

pattern_list <- c(7,8)
size_list <- c(200)


#specific function for estimate T
getT <- function(adj.matrix){
  #for mean degree
  Tdegree = sum(adj.matrix)/NROW(adj.matrix)
  
  #for clustering coefficient 
  Tcc = (transitivity(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))
  
  #for mean betweenness
  Tmb = mean(betweenness(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))
  
  #tg
  T_tg = sum(count_triangles(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))/(choose(NROW(adj.matrix),3)*3)
  return(c(Tdegree,Tcc,Tmb,T_tg))
}

#tool function for parallel
generate_graph <- function(x,size){
  if(x<=2){
    #graph1-2
    W <- graphon_u_to_p_smoothness2(size=size,smoothness = c(0.01,10)[x],sampling_on_u = node_sampling)
  }else if(x<=4){
    #graph3-4
    #W <- graphon_u_to_p(size=size,pattern = "SAS6",smoothness = c(1000,1.5)[x-2],sampling_on_u = node_sampling)
    W <- graphon_u_to_p_smoothness(size=size,smoothness = c(8,2)[x-2])
  }else if(x<=6){
    #graph5-6
    W <- graphon_u_to_p(size=size,pattern = "sinsin",smoothness= c(8,10)[x-4],sampling_on_u = node_sampling)
  }else if(x == 7){
    #read real data
    rat_graph <- igraph::read_graph("../read_data/rattus.norvegicus_brain_3.graphml",format ="graphml")
    rat_graph <- igraph::as.undirected(rat_graph, mode =  "each")
    pop.adj <- as.matrix(igraph::get.adjacency(rat_graph))
    W <- ifelse(pop.adj==2,1,pop.adj)
  }else if(x == 8){
    eu_email_graph <- igraph::read_graph("../read_data/email-Eu-core.txt")
    eu_email_graph <- igraph::as.undirected(eu_email_graph, mode =  "each")
    pop.adj <- as.matrix(igraph::get.adjacency(eu_email_graph))
    W <- ifelse(pop.adj==2,1,pop.adj)
  }
  
  return(W)
}

library(blockmodels)
#compute se
time_start <- Sys.time()
EstT_list <- mclapply(1:M, function(m) {
  library(tictoc)
  tic()
  se_est_array <- array(NA,c(length(pattern_list),length(size_list),4))
  for(pattern in pattern_list){
    for(size in size_list){
      if(pattern == 3){
        best_nb = c(0.4,0.1,0.1,0.1)[which(size_list==size)]
      }else{
        best_nb = c(0.1,0.1,0.1,0.1)[which(size_list==size)]
      }
      
      
      W <- generate_graph(pattern,size)
      adj.matrix <- gmodel.P(W,symmetric.out = TRUE)
      
      #blist <- sample(1:NROW(W),size,replace = TRUE)
      #adj.matrix <- W[blist,blist]
      
      #library(tictoc)
      #tic()
      my_model <- BM_bernoulli("SBM",adj.matrix,verbosity=0,plotting = "")
      my_model$estimate()
      #toc()
      best_nblock = which.max(my_model$ICL)
      i_model = my_model$memberships[[best_nblock]]
      Z = i_model$Z
      memberships = apply(Z,1,which.max)
      p = my_model$model_parameters[[best_nblock]]$pi
      est_p = p[memberships,memberships]
      Tresult = matrix(NA,B,4)
      for(b in 1:B){
        boot_adj = P_nb_boot(est_p,1)[[1]]
        Tresult[b,] = getT(boot_adj)
      }
      Tresult <- apply(Tresult,2,sd)
      se_est_array[which(pattern_list==pattern),which(size_list==size),] <- Tresult
    }
  }
  toc()
  return(se_est_array)
}, mc.cores=cl_nodes) 
time_end <- Sys.time()
print(time_end-time_start)

result_array_boot_all = array(NA,c(length(pattern_list),length(size_list),4,M))
m=1
for(i in (1:length(EstT_list))){
  result_array_boot_all[,,,m] = EstT_list[[i]]
  m=m+1
}
saveRDS(result_array_boot_all,"results_hist_d1graph_d2_size_d3type_d4sim_realdatasim.RDS")
true_stats_boot = apply(result_array_boot_all,c(1,2,3),mean)
true_stats_boot

size = 800
M = 20
#M_true = 2400
B = 500
node_sampling = TRUE
cl_nodes = 10

pattern_list <- c(3,5)
size_list <- c(200,400,800,2000)


#specific function for estimate T
getT <- function(adj.matrix){
  #for mean degree
  #sum(adj.matrix)/NROW(adj.matrix)
  
  #for clustering coefficient 
  #(transitivity(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))
  
  #for mean betweenness
  #mean(betweenness(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))
  
  #tg
  sum(count_triangles(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))/(choose(NROW(adj.matrix),3)*3)
  
}

#tool function for parallel
generate_graph <- function(x,size){
  if(x<=2){
    #graph1-2
    W <- graphon_u_to_p_smoothness2(size=size,smoothness = c(0.01,10)[x],sampling_on_u = node_sampling)
  }else if(x<=4){
    #graph3-4
    #W <- graphon_u_to_p(size=size,pattern = "SAS6",smoothness = c(1000,1.5)[x-2],sampling_on_u = node_sampling)
    W <- graphon_u_to_p_smoothness(size=size,smoothness = c(8,2)[x-2])
  }else if(x<=6){
    #graph5-6
    W <- graphon_u_to_p(size=size,pattern = "sinsin",smoothness= c(8,10)[x-4],sampling_on_u = node_sampling)
  }
  return(W)
}


#generate true se
time_start <- Sys.time()
TrueT_list <- mclapply(1:M, function(m) {
  se_est_v <- c() 
  for(pattern in pattern_list){
    for(size in size_list){
      W <- generate_graph(pattern,size)
      adj.matrix <- gmodel.P(W,symmetric.out = TRUE)
      library(tictoc)
      #tic()
      boot.Tlist <- zhu_nb_boot(adj.matrix,0.1,B,kowning_u=NULL,method="own",
                                induced_sampling = TRUE,weighted=FALSE,fast=1,getT = getT)
#toc()

se_est <- sd(unlist(boot.Tlist))

T <- getT(adj.matrix)
se_est_v <- c(se_est_v,se_est)
  }
}
  return(se_est_v)
  }, mc.cores=cl_nodes) 
time_end <- Sys.time()
print(time_end-time_start)
apply(do.call("cbind",TrueT_list),1,mean)





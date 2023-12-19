# -------------------------------------------------------------------
# Simulation Script: True Standard Errors for Graph Statistics
# -------------------------------------------------------------------
# Description:
# This script is designed to generate simulated various network and obtain
# true standard errors for various graph statistics, including global 
# clustering coefficient, triangle density, mean degree, and mean betweenness.
# It covers various types of networks generated under different settings.
#
# This script can be used to produce in Table 1 of xxxxx. 
#
# The script will:
# 1. Generate networks using predefined parameters and patterns.
# 2. Compute graph statistics for each network.
# 3. Estimate the true standard errors for these statistics using
#    a very high number of simulations.
#
# Networks include both synthetic and mimicking real-world data examples.
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
# Script Starts Here
# -------------------------------------------------------------------


library(parallel)
library(foreach)
library(igraph)
#library(viridis)
source('../graph_utils.r')
source('../graph_boot_funs.R')
#Sys.sleep(3600*5)
size = 200
#M = 100
M_true = 100000
B = 500
node_sampling = TRUE
cl_nodes = 12
#run for sim real data for now
pattern_list <- c(7,8)#c(1,2,3,4,5,6)
size_list <- c(200)


#specific function for estimate T
getT <- function(adj.matrix){
  #for mean degree
  #sum(adj.matrix)/NROW(adj.matrix)
  
  #for clustering coefficient 
  #(transitivity(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))
  
  #for mean betweenness
  #mean(betweenness(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))
  
  #tg
  #sum(count_triangles(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))/sum(count_triangles(make_full_graph(NROW(adj.matrix))))
  
  #all four
  c(sum(adj.matrix)/NROW(adj.matrix),
    (transitivity(graph_from_adjacency_matrix(adj.matrix,mode = "undirected"))),
    mean(betweenness(graph_from_adjacency_matrix(adj.matrix,mode = "undirected"))),
    sum(count_triangles(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))/sum(count_triangles(make_full_graph(NROW(adj.matrix))))
  )
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


#generate true se
time_start <- Sys.time()
TrueT_list <- mclapply(1:M_true, function(m) {
  true_T <- matrix(NA,0,4) 
  for(pattern in pattern_list){
    #for(size in size_list){
    W <- generate_graph(pattern,size)
    
    if(pattern %in% c(7,8)){
      blist <- sample(1:NROW(W),size,replace = TRUE)
      adj.matrix <- W[blist,blist]
    }else{
      adj.matrix <- gmodel.P(W,symmetric.out = TRUE)
    }
    
    T <- getT(adj.matrix)
    true_T <- rbind(true_T,T)
    #}
  }
  return(true_T)
}, mc.cores=cl_nodes) 
time_end <- Sys.time()
print(time_end-time_start)
result_array = array(unlist(TrueT_list),dim=c(length(pattern_list),4,M_true))
apply(result_array,c(1,2),sd)

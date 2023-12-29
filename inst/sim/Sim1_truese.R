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
# Created on: 12/28/2023
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
library(localboot) 

# -------------------------------------------------------------------
# Script Configuration
# -------------------------------------------------------------------

# size of the network, number of node
size = 200

# number of simulations 
M_true = 10000 # use number greater than 500,000 for real simulation

# sampling on node indices
node_sampling = TRUE

# define a tool function to obtain graph statistics of interest
# for example, we start with clustering coefficient
getT <- function(adj.matrix){
  #for clustering coefficient 
  (transitivity(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))
}

# For other graph statistics:
#getT <- function(adj.matrix){
  #for mean degree
  #sum(adj.matrix)/NROW(adj.matrix)
  #for clustering coefficient 
  #(transitivity(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))
  #for mean betweenness
  #mean(betweenness(graph_from_adjacency_matrix(adj.matrix,mode = "undirected"))
  #triangle density
  #sum(count_triangles(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")))/sum(count_triangles(make_full_graph(NROW(adj.matrix))))
#}

# For parallel computing, number of CPU cores
cl_nodes = 8

# -------------------------------------------------------------------
# Script Starts Here
# -------------------------------------------------------------------

#specify network type
pattern_list <- c(1,2,3,4,5,6) # c(7,8) for real data

#generate true se
time_start <- Sys.time()
TrueT_list <- mclapply(1:M_true, function(m) {
  true_T <- c()
  for(pattern in pattern_list){
    P <- generate_graphon(size,pattern)
    adj.matrix <- localboot::generate_network_P(P)
    Graph_T <- getT(adj.matrix)
    true_T <- c(true_T,Graph_T)
  }
  return(true_T)
}, mc.cores=cl_nodes) 
time_end <- Sys.time()
print(time_end-time_start)
result_array = array(unlist(TrueT_list),dim=c(length(pattern_list),M_true))
apply(result_array,1,sd)

# -------------------------------------------------------------------
# Simulation Script: Results from Table 1, Estimated Standard Errors for Graph Statistics
# -------------------------------------------------------------------
# Description:
# This script is designed to generate simulated various network and obtain estimated
# standard errors for various graph statistics, including global 
# clustering coefficient, triangle density, mean degree, and mean betweenness.
# It covers various types of networks generated under different settings.
# 
# This script can be used to produce in Table 1 of xxxxx. 
#
#
# The script will:
# 1. Set up the necessary parameters as outlined in Table 1.
# 2. Apply various methods including proposed methods on various 
# generated networks.
# 3. Report the average estimated standard errors for selected graph statistics. 
#
# Note:
# Due to the complexity and size of the data in Table 1, running this
# script may require significant computational resources and time with the actual
# replication number and network size reported in Table 1.
# Users should be prepared for long execution times and ensure
# sufficient system resources are available.
#
# Author: Tianhai Zu
# Affiliation: University of Texas at San Antonio
# Created on: 12/18/2023
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
M = 50 # use number greater than 1,000 for real simulation

# number of bootstrap networks generated each time
B = 500

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
TrueT_list <- mclapply(1:M, function(m) {
  true_T <- c()
  for(pattern in pattern_list){
    P <- generate_graphon(size,pattern)
    adj.matrix <- localboot::generate_network_P(P)
    #fit local boot
    local_boot_res = local_boot(adj.matrix,B,returns = "T",getT=getT)
    Graph_T <- sd(unlist(local_boot_res))
    true_T <- c(true_T,Graph_T)
  }
  return(true_T)
}, mc.cores=cl_nodes) 
time_end <- Sys.time()
print(time_end-time_start)
result_array = array(unlist(TrueT_list),dim=c(length(pattern_list),M_true))
apply(result_array,1,mean)

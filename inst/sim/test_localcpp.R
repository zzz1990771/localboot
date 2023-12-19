# -------------------------------------------------------------------
# Internal Script: Test CPP accelerated localboot
# -------------------------------------------------------------------

#
# Author: Tianhai Zu
# Affiliation: University of Texas at San Antonio
# Created on: 12/13/2023
#
# Requirements:
# The script requires the following R packages: 
# - igraph (for network generation and manipulation)
# - localboot (this package)
# - ggplot2 (for data visualization)
# Install these packages using install.packages() if not already installed.

# Library Imports
library(igraph)
#library(localboot) #currently commented out
library(ggplot2)

#temp source, will be deleted after the package is finished.
source("./R/generate_network.R")
source("./R/localboot.R")
source("./R/plot_utils.R")

# -------------------------------------------------------------------
# Script Starts Here
# -------------------------------------------------------------------

# Number of node in a network
n = 400

# Generate network adjacency matrix
P <- generate_graphon(size = n, graph_num = 1)
adj.matrix <- generate_network_P(P)

# Let's assume we are investigating:
# **global clustering coefficient**

# We define a tool function to specify the target graph statistics 
getT <- function(adj.matrix){
  library(igraph)
  #clustering coefficient 
  return( transitivity(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")) )
}

# obtain global clustering coefficient from the generated network
getT(adj.matrix)

# should be similar
localboot    (adj.matrix,0.1,200,returns = "T")$se 
localboot_old(adj.matrix,0.1,200,returns = "T")$se

# benchmark run time without returning graph statistics
library(microbenchmark)
benchmark_result <- microbenchmark(
  function1_run = localboot(adj.matrix,0.1,200),  
  function2_run = localboot_old(adj.matrix,0.1,200),  
  times = 10
)
benchmark_summary <- summary(benchmark_result)
print(benchmark_summary)

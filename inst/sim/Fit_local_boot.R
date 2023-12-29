# -------------------------------------------------------------------
# Demonstration Script: Local Bootstrap on Generated Network
# -------------------------------------------------------------------
# Description:
# This script demonstrates how to apply a local bootstrap procedure
# to a network generated with our functions. It is designed to showcase
# methods for bootstrapping in network analysis, particularly
# obtaining standard error estimates.
#
# The script will:
# 1. Generate a network using predefined parameters.
# 2. Implement the local bootstrap procedure on this network.
# 3. Report the estiamted standard error on global clustering coefficient.
#
#
# Author: Tianhai Zu
# Affiliation: University of Texas at San Antonio
# Created on: 12/05/2023
#
# Requirements:
# The script requires the following R packages: 
# - igraph (for network generation and manipulation)
# - localboot (this package)
# - ggplot2 (for data visualization)
# Install these packages using install.packages() if not already installed.

# Library Imports
library(igraph)
library(localboot)
library(ggplot2)

# -------------------------------------------------------------------
# Script Starts Here
# -------------------------------------------------------------------

# Number of node in a network
n = 400

# Generate network adjacency matrix
P <- generate_graphon(n,1) 
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

# Use local bootstrap to estimate standard error of certain graph statistics
local_boot_res = localboot(adj.matrix,200,returns = "T",getT=getT)
local_boot_res$se #should be around 0.003 to 0.0036.

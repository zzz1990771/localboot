# -------------------------------------------------------------------
# Demostration Script: Graphical Models and Graphons
# -------------------------------------------------------------------
# Description:
# This script is designed to generate a variety of networks using
# different size and network settings. It serves as a utility for creating 
# synthetic networks to test network analysis algorithms.
#
# The script will:
# 1. Generate networks based on predefined graphical models.
# 2. Visualize the graphon probability matrix and adjacency matrix.
# 3. Use getT() function to obtain various estimated graph statistics. 
# 4. Save the generated networks for further analysis or use in
#    simulation studies.
#
# This script can be modified to accommodate different network
# sizes and different network settings as needed for specific research
# or teaching purposes. We provide six network settings with the package locaboot.
#
# Author: Tianhai Zu
# Affiliation: University of Texas at San Antonio
# Created on: 12/18/2023
#
# Requirements:
# The script requires the following R packages: 
# - igraph (for network generation and analysis)
# - localboot (this package)
# Install these packages using install.packages() if not already installed.
# Additional files, like real-world network data, should be in the
# specified directory.

# Library Imports
library(igraph)
#library(localboot) #currently commented out

# temp source, will be deleted
source('./R/generate_network.R')

# -------------------------------------------------------------------
# Script Configuration
# -------------------------------------------------------------------

# Number of node in a network
n = 400

# Specify network type
# There are six network settings labeled from 1 - 6
# Please feel free to try each setting and use the following code to
# visualize the probability matrix and adjacency matrix, and use getT()
# to obtain various estimated graph statistics on this network.
graphon_no = 1

# To generate simulated network from real network, specify graphon_no as 7 or 8.
# 7:
# 8:
# graphon_no = 7


# -------------------------------------------------------------------
# Script Starts Here
# -------------------------------------------------------------------

# Generate graphon probability matrix 
P <- generate_graphon(size = n, graph_num = graphon_no)

# Generate network adjacency matrix
adj.matrix <- generate_network_P(P)

# Visualize the graphon P matrix
plot_P(P)
# Visualize the graphon adjacency matrix
plot_adj(adj.matrix)

# Estimate global clustering coefficient
getT <- function(adj.matrix){
  library(igraph)
  #clustering coefficient 
  return( transitivity(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")) )
}
getT(adj.matrix)
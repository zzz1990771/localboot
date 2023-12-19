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
#library(localboot) #currently commented out

# temp source, will be deleted
source('./R/')
source('../graph_boot_funs.R')

# -------------------------------------------------------------------
# Script Configuration
# -------------------------------------------------------------------

# Number of node in a network
n = 400

# Specify network type
graphon_no = 1

# Let's assume we are investigating:
# **global clustering coefficient**
# We define a tool function to specify the target graph statistics 
getT <- function(adj.matrix){
  library(igraph)
  #clustering coefficient 
  return( transitivity(graph_from_adjacency_matrix(adj.matrix,mode = "undirected")) )
}

# -------------------------------------------------------------------
# Script Starts Here
# -------------------------------------------------------------------


# Generate network adjacency matrix
P <- generate_graphon(n,graphon_no)
adj.matrix <- generate_network_P(P)




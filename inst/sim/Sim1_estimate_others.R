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
# C. Embedding Method - Employs embedding to bootstrap networks.
# G. Histogram Bootstrap - Implements a histogram-based bootstrap approach for network data.
# 
# Methods not provided here as they do not really work:
# D. Neighborhood-Sampling Method
# E. Proposed Method with Edge Sampling
# F. Proposed Method with Random Neighbors
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
# - blockmodels (package for histogram method)
# - RSpectra (package for subsampling method)
# Install these packages using install.packages() if not already installed.
# Additional files, like real-world network data, should be in the
# specified directory.

# Library Imports
library(igraph)
library(parallel)
library(foreach)
library(localboot) 

#for other methods in this simulation script, we need to load
source(system.file("sim", "sup_funcs_sim.R", package = "localboot"))
library(blockmodels)
library(RSpectra)

# -------------------------------------------------------------------
# Script Configuration
# -------------------------------------------------------------------

# size of the network, number of node
size = 800

# number of simulations 
M = 18*3 # use number greater than 1,000 for real simulation

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
cl_nodes = 18

#specify network type
pattern_list <- c(1,2,3,4,5,6) # c(7,8) for real data
pattern_list <- c(1,2) # c(7,8) for real data

# -------------------------------------------------------------------
# Script Starts Here
# -------------------------------------------------------------------

# Method A: Empirical Graphon Bootstrap (EGB)
# As this method is a special case of local bootstrap, we use 
# function `localboot` with size of neighbors being 0 to achieve EGB.

# obtain estimated standard errors via EGB
time_start <- Sys.time()
EGB_list <- mclapply(1:M, function(m) {
  est_T <- c()
  for(pattern in pattern_list){
    P <- generate_graphon(size,pattern)
    adj.matrix <- localboot::generate_network_P(P)
    #fit EGB via localboot
    EGB_boot_res = localboot(adj.matrix,B,quantile_n = 1/size, returns = "T",getT=getT)
    Graph_T <- EGB_boot_res$se
    est_T <- c(est_T,Graph_T)
  }
  return(est_T)
}, mc.cores=cl_nodes) 
time_end <- Sys.time()
print(time_end-time_start)
result_array = array(unlist(EGB_list),dim=c(length(pattern_list),M))

#print out
print("Method A: Empirical Graphon Bootstrap (EGB)")
print("Estimated SE for various networks:")
print(apply(result_array,1,mean))

# Method B: subsampling 

for(sample_size in seq(0.95,0.99,0.01)){
  # obtain estimated standard errors via subsampling
  time_start <- Sys.time()
  sub_list <- mclapply(1:M, function(m) {
    est_T <- c()
    for(pattern in pattern_list){
      P <- generate_graphon(size,pattern)
      adj.matrix <- localboot::generate_network_P(P)
      #fit subsampling boot
      ls.boot.num <- vertex_subsampling(adj.matrix,b=sample_size*size,N=B,my_function=getT,output_size = 1)
      ls.se <- sqrt(subsampling_var.numeric(numeric=ls.boot.num,n =size,b=sample_size*size))
      Graph_T <- ls.se
      est_T <- c(est_T,Graph_T)
    }
    return(est_T)
  }, mc.cores=cl_nodes)
  time_end <- Sys.time()
  print(time_end-time_start)
  result_array = array(unlist(sub_list),dim=c(length(pattern_list),M))
  #print out
  print("Method B: subsampling")
  print("Estimated SE for various networks:")
  print(sample_size)
  print(apply(result_array,1,mean))
}

# obtain estimated standard errors via subsampling
time_start <- Sys.time()
sub_list <- mclapply(1:M, function(m) {
  est_T <- c()
  for(pattern in pattern_list){
    P <- generate_graphon(size,pattern)
    adj.matrix <- localboot::generate_network_P(P)
    #fit subsampling boot
    ls.boot.num <- vertex_subsampling(adj.matrix,b=0.1*size,N=500,my_function=getT,output_size = 1)
    ls.se <- sqrt(subsampling_var.numeric(numeric=ls.boot.num,n=size,b=0.1*size))
    Graph_T <- ls.se
    est_T <- c(est_T,Graph_T)
  }
  return(est_T)
}, mc.cores=cl_nodes)
time_end <- Sys.time()
print(time_end-time_start)
result_array = array(unlist(sub_list),dim=c(length(pattern_list),M))
#print out
print("Method B: subsampling")
print("Estimated SE for various networks:")
apply(result_array,1,mean)



# Method C: Embedding
# obtain estimated standard errors via Embedding
time_start <- Sys.time()
EstT_list <- mclapply(1:M, function(m) {
  est_T <- c()
  for(pattern in pattern_list){
    P <- generate_graphon(size,pattern)
    adj.matrix <- localboot::generate_network_P(P)
    #fit Embedding boot
    local_boot_res = ll_boot(adj.matrix,B,d=5)
    Graph_T <- sd(unlist(lapply(local_boot_res,getT)))
    est_T <- c(est_T,Graph_T)
  }
  return(est_T)
}, mc.cores=cl_nodes) 
time_end <- Sys.time()
print(time_end-time_start)
result_array = array(unlist(EstT_list),dim=c(length(pattern_list),M))

#print out
print("Method C: Embedding")
print("Estimated SE for various networks:")
apply(result_array,1,mean)

# Method: Histogram Bootstrap

# obtain estimated standard errors via Histogram Bootstrap
time_start <- Sys.time()
EstT_list <- mclapply(1:M, function(m) {
  est_T <- c()
  for(pattern in pattern_list){
    P <- generate_graphon(size,pattern)
    adj.matrix <- localboot::generate_network_P(P)
    #step1: get Histogram est_p
    my_model <- BM_bernoulli("SBM",adj.matrix,verbosity=0,plotting = "")
    my_model$estimate()
    best_nblock = which.max(my_model$ICL)
    i_model = my_model$memberships[[best_nblock]]
    Z = i_model$Z
    memberships = apply(Z,1,which.max)
    p = my_model$model_parameters[[best_nblock]]$pi
    est_p = p[memberships,memberships]
    #step2: generate Histogram bootstrap networks
    boot_networks = NULL
    for(b in 1:B){
      resample_index = sample(1:size,size,replace = T)
      resampled_est_p = est_p[resample_index,resample_index]
      boot_networks[[b]] = generate_network_P(resampled_est_p)
    }
    Graph_T <- sd(unlist(lapply(boot_networks,getT)))
    est_T <- c(est_T,Graph_T)
  }
  return(est_T)
}, mc.cores=cl_nodes) 
time_end <- Sys.time()
print(time_end-time_start)
result_array = array(unlist(EstT_list),dim=c(length(pattern_list),M))

#print out
print("Method G: Histogram Bootstrap")
print("Estimated SE for various networks:")
apply(result_array,1,mean)



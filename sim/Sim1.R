# -------------------------------------------------------------------
# Simulation Script: Results from Table 1
# -------------------------------------------------------------------
# Description:
# This script is designed to replicate the results presented in Table 1
# of manuscript xxxxxxxxxx. Table 1 contains extensive data
# and analyses, making this replication process intricate and
# potentially time-consuming.
#
# The script will:
# 1. Set up the necessary parameters as outlined in Table 1.
# 2. Apply various methods including proposed methods on various 
# generated networks.
# 3. Display the results similar to Table 1 but less replication and 
# lower network size.
#
# Note:
# Due to the complexity and size of the data in Table 1, running this
# script may require significant computational resources and time.
# Users should be prepared for long execution times and ensure
# sufficient system resources are available.
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


# -------------------------------------------------------------------
# Script Starts Here
# -------------------------------------------------------------------

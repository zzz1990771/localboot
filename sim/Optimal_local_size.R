# -------------------------------------------------------------------
# Plot Generation Script: Optimal Neighbor Set Size for SBMs
# -------------------------------------------------------------------
# Description:
# This script generates plots to visualize the optimal neighbor set size
# for standard error estimation in relation to the number of blocks (K)
# in a Stochastic Block Model (SBM). The goal is to identify and
# illustrate how the choice of neighbor set size impacts the precision
# of standard error estimates in different block model configurations.
#
# The script will:
# 1. Simulate SBM networks with varying numbers of blocks (K).
# 2. Estimate standard errors for each configuration using different
#    neighbor set sizes.
# 3. Generate plots showing the relationship between the optimal neighbor
#    set size and the number of blocks in the SBM.
#
# These visualizations aid in understanding the dynamics of neighbor
# selection in block models and guide the selection of appropriate
# parameters for accurate error estimation.
#
# Author: Tianhai Zu
# Affiliation: University of Texas at San Antonio
# Created on: 12/05/2023
#
# Requirements:
# The script requires the following R packages:
# - igraph (for basic network generation and manipulation)
# - [Any other packages needed for specific models or patterns]
# Install these packages using install.packages() if not already installed.

# Library Imports
library(igraph)
# Additional library imports as required

# Library Imports
library(igraph)
library(ggplot2)
# Additional library imports as required

# -------------------------------------------------------------------
# Script Configuration
# -------------------------------------------------------------------
# Set parameters, such as ranges for K, neighbor set sizes,
# and other simulation settings here.



# -------------------------------------------------------------------
# Script Starts Here
# -------------------------------------------------------------------

# plot-estimation-vs-bootstrap --------------------------------------------
library(ggplot2)
library(readr)
er_g4_g10_cc <- read_csv("er_g4_g10_cc.csv")
#View(er_g4_g10_cc)

#p
plot_er <- ggplot(er_g4_g10_cc, aes(x=grid)) +
  geom_line( aes(y=abs_diff_se_er), size=1, color="blue") + 
  geom_line( aes(y=(f_diff_P_er/600000)), size=1, color="red") +
  scale_y_continuous(
    # Features of the first axis frac(1,m) 
    name = expression( abs(widehat(S.E.) - S.E.)),#,
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*600000, name="",)
  ) + scale_x_continuous(name = "s/n") + 
  theme(text = element_text(size=18),
        axis.text.y.right = element_text(colour="red"), 
        axis.text.y.left = element_text(colour="blue"))
#theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
#      axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
#      axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
#      axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"))

plot_g4 <- ggplot(er_g4_g10_cc, aes(x=grid)) +
  geom_line( aes(y=abs_diff_se_g4), size=1, color="blue") + 
  geom_line( aes(y=(f_diff_P_g4/200000)), size=1, color="red") +
  scale_y_continuous(
    # Features of the first axis
    name = " ",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*200000, name=" "),
    labels = scales::scientific
  ) + scale_x_continuous(name = "s/n") + 
  theme(text = element_text(size=18),
        axis.text.y.right = element_text(colour="red"), 
        axis.text.y.left = element_text(colour="blue"))

plot_g10 <- ggplot(er_g4_g10_cc, aes(x=grid)) +
  geom_line( aes(y=abs_diff_se_g10), size=1, color="blue") + 
  geom_line( aes(y=(f_diff_P_10/600000)), size=1, color="red") +
  scale_y_continuous(
    # Features of the first axis
    name = " ",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*600000, name=expression(abs(widehat(P) - P)[2]^2))
  ) + scale_x_continuous(name = "s/n") + 
  theme(text = element_text(size=18),
        axis.text.y.right = element_text(colour="red"), 
        axis.text.y.left = element_text(colour="blue"))

require(gridExtra)
grid.arrange(plot_er, plot_g4, plot_g10,ncol=3)

# plot-bestnv-vs-K with ev=0 or not --------------------------------------------
library(ggplot2)
library(tidyverse)
library(readr)
nb_vs_K <- read_csv("nb_vs_K.csv")

#oracle
nb_vs_K_oracle <- filter(nb_vs_K,dist=="oracle")
nb_vs_K_oracle$best_nb_min_sma3 <- NA
nb_vs_K_oracle$best_nb_min_sma3[1:19] <- TTR::SMA(nb_vs_K_oracle$best_nb_min[1:19],4)
nb_vs_K_oracle$best_nb_min_sma3[1:3] <- nb_vs_K_oracle$best_nb_min[1:3]

nb_vs_K_oracle$best_nb_min_sma3[20:38] <- TTR::SMA(nb_vs_K_oracle$best_nb_min[20:38],4)
nb_vs_K_oracle$best_nb_min_sma3[20:22] <- c(0.65,0.5,0.4)

plot_oracle <- ggplot(nb_vs_K_oracle, aes(x=log2(K),group=p,colour=p)) + 
  geom_line( aes(y=best_nb_min_sma3),size=1.2) + 
  scale_y_continuous(name = " best s/n") + 
  #geom_smooth( aes(y=best_nb_min), size=1,se = FALSE) + 
  scale_x_continuous(n.breaks = 12, labels = function(x) { 2^x}, name = "K") +
  theme(legend.title=element_blank(),text = element_text(size=20))

#up
nb_vs_K_up <- filter(nb_vs_K,dist=="up")
nb_vs_K_up$best_nb_min_sma3 <- NA
nb_vs_K_up$best_nb_min_sma3[1:19] <- TTR::SMA(nb_vs_K_up$best_nb_min[1:19],4)
nb_vs_K_up$best_nb_min_sma3[1:3] <- nb_vs_K_up$best_nb_min[1:3]

nb_vs_K_up$best_nb_min_sma3[20:38] <- TTR::SMA(nb_vs_K_up$best_nb_min[20:38],4)
nb_vs_K_up$best_nb_min_sma3[20:22] <- c(0.7,0.58,0.5)

nb_vs_K_up$best_nb_min_sma3 <- nb_vs_K_up$best_nb_min_sma3*200

plot_up <- ggplot(nb_vs_K_up, aes(x=log2(K),group=p,colour=p)) + 
  geom_line( aes(y=best_nb_min_sma3),size=1.2) + 
  scale_y_continuous(name = "Optimal neighbor set size s") + 
  #geom_smooth( aes(y=best_nb_min), size=1,se = FALSE) + 
  scale_x_continuous(n.breaks = 12, labels = function(x) { 2^x}, name = "K") +
  theme(legend.title=element_blank(),text = element_text(size=20))


require(gridExtra)
grid.arrange(plot_oracle, plot_up,ncol=2)

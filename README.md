# localboot

## Overview

`localboot` is an R package that offers tools for bootstrapping various networks through a local bootstrap procedure. It's particularly useful in network analysis for estimating the uncertainty of graph statistics. The package also includes utility functions that generating probability matrices, creating network adjacency matrices, and plotting network.

## Installation

To install the latest version of `localboot` from CRAN:

```R
install.packages("localboot")
```

To install the development version from GitHub:

```R
# install.packages("devtools")
devtools::install_github("zzz1990771/localboot")
```

## Features

- Local bootstrap methods for network analysis.
- Tools for generating probability matrices and network adjacency matrices.
- Functions for plotting network structures.
- Several simulation scripts included in the `inst\sim` folder for advanced analysis.

## Usage

After installing the `localboot` package, you can load it into your R session:

```R
library(localboot)
?localboot
```

Example usage:

```R
# Example usage
P = generate_graphon(100, 1)
A = generate_network_P(P, replicate = 1, symmetric.out = TRUE)
result <- localboot(A = A, B = 100, returns = "boot")
```

### Included Simulation Scripts

The localboot package comes with several simulation scripts located in the `inst\sim` folder. These scripts can be loaded and used for various simulations and analyses:

- `Fit_local_boot.R`: Demonstrates the application of local bootstrap procedures to networks. It generates random networks by `generate_graphon()` and `generate_network_P()` and fit local boostrap procedures on it using the main function `localboot()`.

- `Generate_various_networks.R`: Generates a variety of synthetic networks for testing network analysis algorithms.

- `Sim1_truese.R`: Obtains true simulated standard errors for graph statistics in various network simulations.

- `Sim1_estimates_local.R`: Generates various simulated networks and obtains estimated standard errors for different graph statistics by applying local bootstrap procedure.

- `Sim1_estimate_others.R`: Generates various simulated networks and obtains estimated standard errors for different graph statistics by applying other alternative procedures.

- `sup_funcs_sim.R`: Support functions for simulations in `Sim1_estimate_others.R`.

- `Optimal_local_size.R`: Generates plots to visualize the optimal neighbor set size of local bootstrap procedure for standard error estimation in relation to the number of blocks (K) in a Stochastic Block Model (SBM). 


To load a script, use the `system.file()` function to find its path and then load it with `file.edit()`. For example:

```R
scriptPath <- system.file("sim/Generate_various_networks.R", package = "localboot")
file.edit(scriptPath)
```

## Results from Paper XXX

This package, along with the scripts `Sim1_truese.R`, `Sim1_estimates_local.R`, and `Sim1_estimate_others.R`, can be used to replicate the results presented in Paper XXX. Detailed instructions and data required for replication are provided within each script, ensuring a comprehensive understanding of the proposed local boostrap procedure.


## License

This package is free and open source software, licensed under GPL-3.

## Authors

-Tianhai Zu (zuti@mail.uc.edu)

-Yichen Qin (qinyn@ucmail.uc.edu)

## LL method, Levin and Levina 2019
ll_boot <- function(A,
                    B,
                    d = 2,
                    induced_sampling = TRUE,
                    getT = NULL,
                    ...) {
  #first ASE method
  #input A adj
  #input d top d largest engenvalues
  ASE <- function(X, d = 2) {
    L = svd(X)
    ud = L$u[, 1:d]
    dd = diag(L$d[1:d])
    ASE_matrix <- ud %*% sqrt(dd)
    ASE_matrix
  }
  
  #basic info
  N=NROW(A)
  
  #replacing A diagonal to scaled degree
  diag(A) <- apply(A, 2, mean)
  
  #obtain X from ASE of A
  X_hat <- ASE(A, d = d)
  
  #empirical distribution of A
  P <- X_hat %*% t(X_hat)
  # limit p to 0-1
  P[P > 1] <- 1
  P[P < 0] <- 0
  
  #generate bootstrap
  boot.list <- list()
  for (i in 1:B) {
    if (induced_sampling) {
      blist <- sample((1:N), N, replace = T)
    } else{
      blist <- 1:N
    }
    
    p_hat_matrix_b <- P[blist, blist]
    random.matrix <- matrix(runif(N * N, 0, 1), N, N)
    random.matrix[lower.tri(random.matrix)] <-
      t(random.matrix)[lower.tri(t(random.matrix))]
    g.adj.nb <- 1 * (random.matrix < p_hat_matrix_b)
    if (is.null(getT)) {
      boot.list[[i]] <- g.adj.nb
    } else{
      boot.list[[i]] <- getT(g.adj.nb)
    }
  }
  return(boot.list)
}


# Lunde Sarkar Codes for paper "Subsampling Sparse Graphons Under Minimal Assumptions"
require(RSpectra)
### Vertex subsampling
# b: subsample size (integer)
# N number of subsamples
vertex_subsampling = function(Adj_mat, b,N, my_function,output_size,...) {
  
  boot_out = matrix(NA,nrow=N, ncol=output_size)
  
  n_row = nrow(Adj_mat)
  
  for(ii in 1: N){
    boot_indices= sample(n_row, size = b, replace=FALSE)
    boot_subgraph = Adj_mat[boot_indices,boot_indices]
    boot_out[ii,] = my_function(boot_subgraph,...) 
  }
  return(boot_out)
}
### P-subsampling
# b: subsample size (integer)
# Subsample probability p
# N number of p-subsamples

p_subsampling = function(Adj_mat,p_n,N,  my_function, output_size,...) {
  boot_out = list(matrix(NA,nrow=N,ncol=output_size), numeric(N))
  
  n_row = nrow(Adj_mat)
  
  for(ii in 1: N){
    # Select nodes with probability p_n
    boot_indices= which(rbinom(n_row,1,prob=p_n) == 1)
    
    # Construct induced subgraph
    Adj_mat_boot = Adj_mat[boot_indices,boot_indices]
    
    # Delete isolated vertices
    included_vertices = which(rowSums(Adj_mat_boot) > 0)
    print(length(included_vertices))
    if(length(included_vertices) >0) {
      boot_subgraph = Adj_mat_boot[included_vertices,included_vertices]
      boot_out[[1]][ii,] = my_function(boot_subgraph,...) 
    } else {
      boot_out[[1]][ii,] = rep(NA,output_size)
    }
    boot_out[[2]][ii] = length(included_vertices)
  }
  return(boot_out)
}

### Functions to compute confidence intervals from subsampling
# est.param: value of size n functional
# value used to center subsample distributions
#  For jackknife bias correction, est.param is set to debiased estimate,
# centering.val is set to subsample means.

subsampling_CI = function(x,...) UseMethod("subsampling_CI")

subsampling_CI.numeric= function(numeric,est.param,centering.val,tau_seq=function(x) return(x^0.5),n,b,alpha=0.05){
  tau_vec = Vectorize(tau_seq,vectorize.args="x")  
  
  quantile_values = c(1-alpha/2,alpha/2)
  subsample_normalized = tau_vec(b)*(numeric-centering.val)
  quantiles_subsample = quantile(subsample_normalized, probs=quantile_values,na.rm=TRUE)
  ci_out = est.param - quantiles_subsample/tau_seq(n)
  return(ci_out)
}

subsampling_CI.matrix= function(matrix,est.param,centering.val,tau_seq=function(x) return(x^0.5),n,b,alpha=0.05){
  num_cols = ncol(matrix)
  num_rows = nrow(matrix)
  
  tau_vec = Vectorize(tau_seq,vectorize.args="x")
  
  quantile_values = c(1-alpha/2,alpha/2)
  centering.mat = matrix(rep(centering.val,num_rows),nrow=num_rows,byrow = TRUE)
  est.param = matrix(rep(est.param,2),nrow=2,byrow = TRUE)
  
  subsample_normalized = tau_vec(b)*(matrix-centering.mat)
  quantiles_subsample = apply(subsample_normalized,2,quantile,probs=quantile_values)
  ci_out = est.param - quantiles_subsample/tau_seq(n) 
  
  return(ci_out)
}

subsampling_var.numeric= function(numeric,tau_seq=function(x) return(x^0.5),n,b,alpha=0.05){
  tau_vec = Vectorize(tau_seq,vectorize.args="x")  
  centering.val = mean(numeric)
  subsample_normalized = tau_vec(b)*(numeric-centering.val)
  subsample_normalized_var = var(subsample_normalized)/tau_seq(n)
  return(subsample_normalized_var)
}

#### Network Jackknife

### Jackknife bias correction
# Takes as input nxn adjacency matrix, an integer specifying length of
# output, a function to apply to the adjacency matrix, and arbitrary inputs
# Returns a list consisting of the jackknife-bias-corrected estimate and
# the original estimate. 

jackknife_bias_correction = function(Adj_mat,output_size,my_function,...){
  n = nrow(adj_mat)
  jack_estimates = matrix(NA, nrow=n,ncol=output_size)
  my_functional = my_function(Adj_mat,...) 
  for(ii in 1:n) {
    print(ii)
    Adj_mat_jack = Adj_mat[-ii,-ii] 
    jack_estimates[ii,] = my_function(Adj_mat_jack,...)
  }
  jack_mean = colMeans(jack_estimates)
  bias_correct = n*my_functional - (n-1)*jack_mean 
  return(list(bias_correct,my_functional))
}

###### Latent Space Bootstrap
# Helper function to simulate from a P matrix fast with Rcpp
#Takes an nxn probability matrix and generates a Bernoulli matrix with
## the corresponding mean
simulate_from_P_mat = function(P_mat,seed=1){
  set.seed(seed)
  cppFunction( 'arma::mat my_function( int n, NumericMatrix P_mat
  ) {
               arma::mat A_mat(n,n,arma::fill::zeros);
               
               for (int ii = 0; ii < n; ii++)
               {
               for(int jj = ii+1; jj < n; jj++) 
               {
               double P_temp = P_mat(ii,jj);
               A_mat(ii,jj) = R::rbinom(1,P_temp);
               }
               }
               A_mat = arma::symmatu(A_mat);
               return(A_mat);
}
',depends = "RcppArmadillo")
  
  output = my_function(n= nrow(P_mat),P_mat=P_mat)
  return(output)
} 


# Helper function to estimate latent positions from observed adjacency matrix.
# Fast eigensolver from RSpectra is used. 
# Inputs: latent_dim: dimension of RDPG model
# eig_method: "which" variable in eigs_sym of RSpectra

extract_latent_positions = function(Adj_mat,latent_dim,eig_method_extract="LM") {
  k_vec = 1:latent_dim
  eigen_decomp = eigs_sym(Adj_mat,k=latent_dim,which=eig_method_extract)
  pos_eigs_index = k_vec[eigen_decomp$values > 0]
  neg_eigs_index = k_vec[-pos_eigs_index]
  
  pos_component = eigen_decomp$vector[,pos_eigs_index]%*% diag(x=eigen_decomp$values[pos_eigs_index]^0.5,nrow= length(pos_eigs_index))
  neg_component = numeric(0)
  
  if(length(neg_eigs_index) == 0) {
    neg_component = matrix(0, nrow=nrow(Adj_mat),ncol=2)
  } else {
    neg_component =eigen_decomp$vector[,neg_eigs_index]%*% diag(x=(-1*eigen_decomp$values[neg_eigs_index])^0.5,nrow= length(neg_eigs_index)) 
  }
  latent_positions = list(pos_component,neg_component)
  return(latent_positions)
}



###### Function to estimate rho
est_rho_function = function(Adj_mat) {
  n=nrow(Adj_mat)
  return(sum(rowSums(Adj_mat))/(n*(n-1)))
}

###### Please read the RSpectra documentation for information about the function inputs
# Function to compute lambda_r(A)/n rho_n.  
# Inputs
# num_eigens: number of eigenvalues to compute
# eig_method:  which variable in eigs_sym from RSpectra.  

eigenvalue_function= function(Adj_mat, num_eigens, eig_method, rho_n){
  
  n_row = nrow(Adj_mat)
  output = eigs_sym(Adj_mat, k = num_eigens, which = eig_method, opts = list(retvec=FALSE,ncv=15))$values
  return(output/(n_row*rho_n))
}
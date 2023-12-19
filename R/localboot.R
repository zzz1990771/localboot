#' Local Bootstrap for Network Data
#'
#' This function applies a local bootstrap method to network data, represented 
#' by an adjacency matrix. It offers various methods and options for bootstrapping,
#' including handling weighted networks and custom distance functions.
#'
#' @param A A square adjacency matrix of the network.
#' @param B The number of bootstrap samples to generate.
#' @param quantile_n The quantile used for neighborhood selection in some methods.
#'     If set to 0 (default), it's calculated as (log(N) / N)^0.5.
#' @param returns Specifies the type of output returned. Possible values are 
#'     "boot" (default), "p_and_time", "p_and_boot", and "T".
#' @param method The method used for bootstrapping. Options are "own" and "zhu".
#' @param dist_func A function to compute the distance matrix. Default is 
#'     `get_dist_default_eigen`.
#' @param kowning_u An optional known 'u' vector for distance calculation.
#' @param induced_sampling A logical indicating whether to use induced sampling.
#'     Defaults to TRUE.
#' @param weighted A logical indicating if the network is weighted. Defaults to FALSE.
#' @param getT An optional function to apply to each bootstrapped sample.
#' @param user_blist An optional user-provided bootstrap list.
#' @param fast A logical indicating if a faster, approximate method should be used. 
#'     Automatically set based on network size if NULL.
#' @param ... Additional arguments passed to other methods.
#'
#' @return Depending on the `returns` argument, this function can return various types 
#'     of outputs including bootstrapped networks, estimated probabilities, computation 
#'     times, and statistics from the `getT` function.
#'
#'
#' @examples
#' # Example usage
#' P = generate_graphon(100, 1)
#' A = generate_network_P(P, replicate = 1, symmetric.out = TRUE)
#' result <- localboot(A = A, B = 100, returns = "boot")
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib localboot
#' @export
localboot <- function(A, B, quantile_n = 0, returns = "boot", method = "own", dist_func = get_dist_default_eigen,
                       kowning_u = NULL, induced_sampling = TRUE, weighted = FALSE, getT = NULL, 
                       user_blist = NULL, fast = NULL, ...) {
  # For timing purpose
  now <- Sys.time()
  
  # Network size
  N <- NROW(A)
  
  # Set 'fast' based on 'N'
  fast <- if(is.null(fast)) N > 400 else fast
  
  # Additional configurations
  no_loop <- sum(diag(A)) == 0
  max_A <- max(A)
  quantile_n <- if(quantile_n == 0) (log(N) / N)^0.5 else quantile_n
  
  #get dist.matrix
  if(is.null(kowning_u)){
    dist.matrix <- dist_func(A)
  }else{
    dist.matrix <- as.matrix(stats::dist(kowning_u))
  }
  
  
  # Main computation block
  if(method == "own") {
    
    #get neighbors
    nb_actual <- ifelse(fast==1,min(30,(ceiling(quantile_n*N))),(ceiling(quantile_n*N))) #(ceiling(quantile_n*N)
    neibors_matrix <- matrix(as.integer(0),nrow = NROW(A),ncol = nb_actual)
    
    for(i in (1:N)){
      neibor_index <- order(dist.matrix[i,], decreasing=FALSE)[1:ceiling(quantile_n*N)]
      if(fast==1){
        sample_fast_index = sample(1:ceiling(quantile_n*N),size = nb_actual)
        neibors_matrix[i,] <- neibor_index[sample_fast_index]
      }else{
        neibors_matrix[i,] <- neibor_index
      }
    }
    #not weighted graph
    if(weighted==FALSE){
      #estimate p for each node pair
      neibors_matrix = neibors_matrix - 1
      mode(neibors_matrix) <- "integer"
      p_hat_matrix <- calculate_p_hat_matrix(A,neibors_matrix)
    }else{
      nb_array <- array(0,c(N,N,max_A+1))
      for(a in 1:(-1+N)){
        for(b in (a):(N)){
          ab_connection <- as.vector(A[neibors_matrix[a,],neibors_matrix[b,]])
          freq_table <- (table(ab_connection)/length(ab_connection))
          full_freq_table <- sapply(0:max_A,function(x){ifelse(sum(names(freq_table)==x)>0,freq_table[names(freq_table)==x],0)})
          nb_array[a,b,] <- cumsum(full_freq_table)
          nb_array[b,a,] <- nb_array[a,b,]
        }
      }
    }
  } else if(method == "zhu") {
    # 3. quantile as logical
    kernel_mat = matrix(0,N,N)
    for (i in 1:N){
      kernel_mat[i,] = as.double(dist.matrix[i,]<stats::quantile(dist.matrix[i,],quantile_n))
    }
    # 4. L1 normalization of each row
    kernel_mat = kernel_mat/(outer(rowSums(kernel_mat),rep(1,N))+1e-10)
    
    # 5. Compute P
    P = kernel_mat %*% A;
    P = (P+t(P))/2;
    p_hat_matrix <- P
  }
  
  finish_time <- Sys.time() - now
  
  # Returning results based on 'returns' parameter
  if(returns == "p_and_time") {
    return(list(p_hat_matrix = p_hat_matrix, finish_time = finish_time))
  }
  
  # Induced sampling to generate bootstrap networks as lists
  nb_boot.list <- list()
  for(i in 1:B) {
    # Sampling computation
    if(induced_sampling){
      blist <- sample((0:(N-1)),N, replace = T) # N-1 as Cpp starts from 0
    }else{
      blist <- 0:(N-1)
    }
    if(!is.null(user_blist)){
      blist <- user_blist
    }
    blist <- sort(blist)
    if(weighted==FALSE){
      rmatrix <- matrix(stats::runif(N*N,0,1),N,N)
      g.adj.nb <- sample_from_p_cpp(p_hat_matrix,blist,rmatrix,no_loop)
      rownames(g.adj.nb) <- colnames(g.adj.nb) <- paste0("v",blist)
      nb_boot.list[[i]]<- g.adj.nb 
    }else{
      nb_array_b <- nb_array[blist,blist,]
      random.matrix <- matrix(stats::runif(N*N,0,1),N,N)
      random.matrix[lower.tri(random.matrix)] <- t(random.matrix)[lower.tri(t(random.matrix))]
      g.adj.nb <- Reduce("+",lapply((1:(max_A+1)),function(x){
        1*(random.matrix>nb_array_b[,,x])
      }))
      nb_boot.list[[i]]<- g.adj.nb 
    }
  }
  
  # Final return based on 'returns' parameter
  if(returns == "p_and_boot") {
    return(list(p_hat_matrix = p_hat_matrix, nb_boot.list = nb_boot.list))
  } else if(returns == "T") {
    boot.Tlist <- lapply(nb_boot.list, getT)
    return(list(boot.Tlist = boot.Tlist, se = stats::sd(unlist(boot.Tlist))))
  } else {
    #return the bootstrap networks as a list
    return(nb_boot.list)
  }
}

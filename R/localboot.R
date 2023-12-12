#main function, local bootstrap
local_boot_old <- function(A,quantile_n=0,B,returns = "boot",method = "own", distance = "zhu",
                        kowning_u=NULL, induced_sampling=TRUE, weighted=FALSE,getT=NULL,user_blist=NULL,fast=NULL,
                        ...){
  
  
  
  now <- Sys.time()
  #SIZE
  N <- NROW(A)
  
  #decide fast level if not specified, for large graph, fast level should be higher
  if(is.null(fast)){
    fast = N > 400
  }
  
  
  #detect for loop
  no_loop <- sum(diag(A)) == 0
  
  #max edge
  max_A <- max(A)
  
  #quantile_n
  if(quantile_n==0){
    quantile_n <- (log(N)/N)^0.5
  }
  
  # dissimilarity measure
  get_dist_zhu <- function(){
    dist.matrix <- matrix(0,nrow = N,ncol = N)
    A_sq = (A%*%A)/N
    for(i in c(1:(N-1))){
      for(ip in c((i+1):N) ){
        tgtvec = abs(A_sq[i,]-A_sq[ip,])
        tgtvec[i] = 0
        tgtvec[ip] = 0
        max_d = max(tgtvec) # tgtvec2 is Li Chen's
        dist.matrix[i,ip] <- max_d
      }
    }
    dist.matrix[lower.tri(dist.matrix)] <- t(dist.matrix)[lower.tri(t(dist.matrix))]
    dist.matrix
  }
  
  get_dist_test <- function(){
    dist.matrix <- matrix(0,nrow = N,ncol = N)
    dist.matrix <- as.matrix(dist(A, method = "euclidean",diag = TRUE))
    dist.matrix
    #dist.matrix[lower.tri(dist.matrix)] <- t(dist.matrix)[lower.tri(t(dist.matrix))]
  }
  if(quantile_n==1){
    p_hat_matrix <- matrix(sum(A)/(N*(N-1)),nrow = NROW(A),ncol = (NROW(A)))
  }else if(method=="own"){
    #get dist
    if(is.null(kowning_u)){
      dist.matrix <- do.call(switch("zhu","zhu"="get_dist_zhu","test"= "get_dist_test"),args = list())
    }else{
      dist.matrix <- as.matrix(dist(kowning_u))
    }
    #get neighbors
    nb_actual <- ifelse(fast==1,min(30,(ceiling(quantile_n*N))),(ceiling(quantile_n*N))) #(ceiling(quantile_n*N)
    neibors_matrix <- matrix(0,nrow = NROW(A),ncol = nb_actual)
    
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
      p_hat_matrix <- matrix(0,nrow = NROW(A),ncol = (NROW(A)))
      for(a in 1:(-1+NROW(p_hat_matrix))){
        for(b in (a+1):(NROW(p_hat_matrix))){
          ab_connection <- A[neibors_matrix[a,],neibors_matrix[b,]]
          p_hat_matrix[a,b] <- mean(ab_connection)
          #p_hat_matrix[a,b] <- sample(as.vector(ab_connection),1)
        }
      }
      for(a in 1:(NROW(p_hat_matrix))){
        b=a
        ab_connection <- A[neibors_matrix[a,],neibors_matrix[b,]]
        #p_hat_matrix[a,b] <- sum(ab_connection)/(N*N-1))
        p_hat_matrix[a,b] <- mean(ab_connection)
        #p_hat_matrix[a,b] <- sample(as.vector(ab_connection),1)
      }
      p_hat_matrix[lower.tri(p_hat_matrix)] <- t(p_hat_matrix)[lower.tri(t(p_hat_matrix))]
    }else{
      #nb_array <- array(0,c(N,N,(ceiling(quantile_n*N)^2)))
      nb_array <- array(0,c(N,N,max_A+1))
      for(a in 1:(-1+N)){
        for(b in (a):(N)){
          ab_connection <- as.vector(A[neibors_matrix[a,],neibors_matrix[b,]])
          freq_table <- (table(ab_connection)/length(ab_connection))
          full_freq_table <- sapply(0:max_A,function(x){ifelse(sum(names(freq_table)==x)>0,freq_table[names(freq_table)==x],0)})
          #print(full_freq_table)
          
          nb_array[a,b,] <- cumsum(full_freq_table)
          nb_array[b,a,] <- nb_array[a,b,]
        }
      }
    }
    
  }else if(method=="zhu"){
    #get dist
    if(is.null(kowning_u)){
      dist.matrix <- do.call(switch("zhu","zhu"="get_dist_zhu","test"= "get_dist_test"),args = list())
    }else{
      dist.matrix <- as.matrix(dist(kowning_u))
    }
    # 3. quantiled as logical
    kernel_mat = matrix(0,N,N)
    for (i in 1:N){
      kernel_mat[i,] = as.double(dist.matrix[i,]<quantile(dist.matrix[i,],quantile_n))
    }
    # 4. L1 normalization of each row
    kernel_mat = kernel_mat/(outer(rowSums(kernel_mat),rep(1,N))+1e-10)
    
    # 5. Compute P
    P = kernel_mat %*% A;
    P = (P+t(P))/2;
    p_hat_matrix <- P
  }
  
  
  finish_time <- Sys.time()-now
  if(returns=="p_and_time"){
    return(list(p_hat_matrix=p_hat_matrix,finish_time=finish_time))
  }
  
  #induced sampling
  nb_boot.list <- list()
  for (i in 1:B){
    if(induced_sampling){
      blist <- sample((1:N),N, replace = T)
    }else{
      blist <- 1:N
    }
    if(!is.null(user_blist)){
      blist <- user_blist
    }
    blist <- sort(blist)
    if(weighted==FALSE){
      p_hat_matrix_b <- p_hat_matrix[blist,blist]
      random.matrix <- matrix(runif(N*N,0,1),N,N)
      random.matrix[lower.tri(random.matrix)] <- t(random.matrix)[lower.tri(t(random.matrix))]
      if(no_loop == TRUE){
        diag(random.matrix) <- 100
      }
      g.adj.nb <- 1*(random.matrix<p_hat_matrix_b)
 
      rownames(g.adj.nb) <- colnames(g.adj.nb) <- paste0("v",blist)
      nb_boot.list[[i]]<- g.adj.nb 

      
    }else{
      nb_array_b <- nb_array[blist,blist,]
      random.matrix <- matrix(runif(N*N,0,1),N,N)
      random.matrix[lower.tri(random.matrix)] <- t(random.matrix)[lower.tri(t(random.matrix))]
      g.adj.nb <- Reduce("+",lapply((1:(max_A+1)),function(x){
        1*(random.matrix>nb_array_b[,,x])
      }))

        nb_boot.list[[i]]<- g.adj.nb 

    }
    
  }
  
  if(returns=="p_and_boot"){
    return(list(p_hat_matrix=p_hat_matrix,nb_boot.list=nb_boot.list))
  }else if(returns == "T"){
    boot.Tlist <- lapply(nb_boot.list,getT)
    return(list(boot.Tlist = boot.Tlist, se = sd(unlist(boot.Tlist))))
  }else{
    return(nb_boot.list)
  }
  
}



#' Local Bootstrap Function
#'
#' @param A A matrix input.
#' @param quantile_n Quantile parameter, default is 0.
#' @param B Integer parameter.
#' @param returns Type of return value, default is "boot".
#' @param method Method to be used, default is "own".
#' @param distance Distance measurement, default is "zhu".
#' @param kowning_u Optional parameter for known user.
#' @param induced_sampling Boolean for induced sampling, default is TRUE.
#' @param weighted Boolean for weighted graph, default is FALSE.
#' @param getT Function to get T, default is NULL.
#' @param user_blist Optional user bootstrap list, default is NULL.
#' @param fast Boolean for fast computation, default is NULL.
#' @param ... Additional arguments.
#' @return Depends on the 'returns' parameter.
#' @export
#'
#' @examples
#' # Example usage
#' result <- local_boot(A = matrix_data, B = 100, returns = "boot")
local_boot <- function(A, quantile_n = 0, B, returns = "boot", method = "own", distance = "zhu",
                       kowning_u = NULL, induced_sampling = TRUE, weighted = FALSE, getT = NULL, 
                       user_blist = NULL, fast = NULL, ...) {
  
  now <- Sys.time()
  N <- NROW(A)
  
  # Set 'fast' based on 'N'
  fast <- if(is.null(fast)) N > 400 else fast
  
  # Additional configurations
  no_loop <- sum(diag(A)) == 0
  max_A <- max(A)
  quantile_n <- if(quantile_n == 0) (log(N) / N)^0.5 else quantile_n
  
  # Internal function to calculate Zhu's distance
  get_dist_zhu <- function() {
    dist.matrix <- matrix(0,nrow = N,ncol = N)
    A_sq = (A%*%A)/N
    for(i in c(1:(N-1))){
      for(ip in c((i+1):N) ){
        tgtvec = abs(A_sq[i,]-A_sq[ip,])
        tgtvec[i] = 0
        tgtvec[ip] = 0
        max_d = max(tgtvec) # tgtvec2 is Li Chen's
        dist.matrix[i,ip] <- max_d
      }
    }
    dist.matrix[lower.tri(dist.matrix)] <- t(dist.matrix)[lower.tri(t(dist.matrix))]
    dist.matrix
  }
  
  # Internal function for test distance
  get_dist_test <- function() {
    dist.matrix <- matrix(0,nrow = N,ncol = N)
    dist.matrix <- as.matrix(dist(A, method = "euclidean",diag = TRUE))
    dist.matrix
  }
  
  # Main computation block
  if(method == "own") {
    #get dist
    if(is.null(kowning_u)){
      dist.matrix <- do.call(switch("zhu","zhu"="get_dist_zhu","test"= "get_dist_test"),args = list())
    }else{
      dist.matrix <- as.matrix(dist(kowning_u))
    }
    #get neighbors
    nb_actual <- ifelse(fast==1,min(30,(ceiling(quantile_n*N))),(ceiling(quantile_n*N))) #(ceiling(quantile_n*N)
    neibors_matrix <- matrix(0,nrow = NROW(A),ncol = nb_actual)
    
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
      p_hat_matrix <- matrix(0,nrow = NROW(A),ncol = (NROW(A)))
      for(a in 1:(-1+NROW(p_hat_matrix))){
        for(b in (a+1):(NROW(p_hat_matrix))){
          ab_connection <- A[neibors_matrix[a,],neibors_matrix[b,]]
          p_hat_matrix[a,b] <- mean(ab_connection)
          #p_hat_matrix[a,b] <- sample(as.vector(ab_connection),1)
        }
      }
      for(a in 1:(NROW(p_hat_matrix))){
        b=a
        ab_connection <- A[neibors_matrix[a,],neibors_matrix[b,]]
        #p_hat_matrix[a,b] <- sum(ab_connection)/(N*N-1))
        p_hat_matrix[a,b] <- mean(ab_connection)
        #p_hat_matrix[a,b] <- sample(as.vector(ab_connection),1)
      }
      p_hat_matrix[lower.tri(p_hat_matrix)] <- t(p_hat_matrix)[lower.tri(t(p_hat_matrix))]
    }else{
      #nb_array <- array(0,c(N,N,(ceiling(quantile_n*N)^2)))
      nb_array <- array(0,c(N,N,max_A+1))
      for(a in 1:(-1+N)){
        for(b in (a):(N)){
          ab_connection <- as.vector(A[neibors_matrix[a,],neibors_matrix[b,]])
          freq_table <- (table(ab_connection)/length(ab_connection))
          full_freq_table <- sapply(0:max_A,function(x){ifelse(sum(names(freq_table)==x)>0,freq_table[names(freq_table)==x],0)})
          #print(full_freq_table)
          
          nb_array[a,b,] <- cumsum(full_freq_table)
          nb_array[b,a,] <- nb_array[a,b,]
        }
      }
    }
  } else if(method == "zhu") {
    #get dist
    if(is.null(kowning_u)){
      dist.matrix <- do.call(switch("zhu","zhu"="get_dist_zhu","test"= "get_dist_test"),args = list())
    }else{
      dist.matrix <- as.matrix(dist(kowning_u))
    }
    # 3. quantiled as logical
    kernel_mat = matrix(0,N,N)
    for (i in 1:N){
      kernel_mat[i,] = as.double(dist.matrix[i,]<quantile(dist.matrix[i,],quantile_n))
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
  
  # Induced sampling
  nb_boot.list <- list()
  for(i in 1:B) {
    # Sampling computation
    if(induced_sampling){
      blist <- sample((1:N),N, replace = T)
    }else{
      blist <- 1:N
    }
    if(!is.null(user_blist)){
      blist <- user_blist
    }
    blist <- sort(blist)
    if(weighted==FALSE){
      p_hat_matrix_b <- p_hat_matrix[blist,blist]
      random.matrix <- matrix(runif(N*N,0,1),N,N)
      random.matrix[lower.tri(random.matrix)] <- t(random.matrix)[lower.tri(t(random.matrix))]
      if(no_loop == TRUE){
        diag(random.matrix) <- 100
      }
      g.adj.nb <- 1*(random.matrix<p_hat_matrix_b)
      
      rownames(g.adj.nb) <- colnames(g.adj.nb) <- paste0("v",blist)
      nb_boot.list[[i]]<- g.adj.nb 
      
      
    }else{
      nb_array_b <- nb_array[blist,blist,]
      random.matrix <- matrix(runif(N*N,0,1),N,N)
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
    return(list(boot.Tlist = boot.Tlist, se = sd(unlist(boot.Tlist))))
  } else {
    return(nb_boot.list)
  }
}

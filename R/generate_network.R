# Define the formula to generate p matrix of graphon1
graphon3 <- function(u, size){
  p_matrix = matrix(0,nrow=size,ncol=size)
  for(i in 1:size){
    for(j in 1:size){
      p <- sin(((u[i]^2+u[j]^2)^(1/4))*0.01)/2 +0.5
      p_matrix[i,j] <- p    
    }
  }
  p_matrix
}

# Define the formula to generate p matrix of graphon2
graphon6 <- function(u, size){
  p_matrix = matrix(0,nrow=size,ncol=size)
  for(i in 1:size){
    for(j in 1:size){
      p <- sin(((u[i]^2+u[j]^2)^(1/4))*10)/2 +0.5
      p_matrix[i,j] <- p    
    }
  }
  p_matrix
}

# Define the formula to generate p matrix of graphon3
graphon1 <- function(u, size){
  p_matrix = matrix(0,nrow=size,ncol=size)
  for(i in 1:size){
    for(j in 1:size){
      p <- cos(pi*(u[i]-u[j]))/(2^3) + 0.5
      p_matrix[i,j] <- p    
    }
  }
  p_matrix
}
# Define the formula to generate p matrix of graphon4
graphon4 <- function(u, size){
  p_matrix = matrix(0,nrow=size,ncol=size)
  for(i in 1:size){
    for(j in 1:size){
      p <- cos(pi*(u[i]-u[j]))/(2) + 0.5
      p_matrix[i,j] <- p    
    }
  }
  p_matrix
}

# Define the formula to generate p matrix of graphon5
graphon2 <- function(u, size){
  p_matrix = matrix(0,nrow=size,ncol=size)
  s = 8
  for(i in 1:size){
    for(j in 1:size){
      if(u[i]>=u[j]){
        p = sin((s*u[i]-s/2)*sin(s*u[j]-s/2))
      }else{
        p = sin((s*u[j]-s/2)*sin(s*u[i]-s/2))
      }
      p_matrix[i,j] <-  (p+1)/2
    }
  }
  p_matrix
}

# Define the formula to generate p matrix of graphon6
graphon5 <- function(u, size){
  p_matrix = matrix(0,nrow=size,ncol=size)
  s = 10
  for(i in 1:size){
    for(j in 1:size){
      if(u[i]>=u[j]){
        p = sin((s*u[i]-s/2)*sin(s*u[j]-s/2))
      }else{
        p = sin((s*u[j]-s/2)*sin(s*u[i]-s/2))
      }
      p_matrix[i,j] <- (p+1)/2
    }
  }
  p_matrix
}

# Function to generate sub-network from real network, labeled as 7
graphon7 <- function(u,size){
  #read real data
  RDS_path <- system.file("rds", "rattus_norvegicus_brain_3_adj.RDS", package = "localboot")
  P <- readRDS(RDS_path)
  node_samples <- sample(1:NROW(P),size,replace = TRUE)
  P[node_samples,node_samples]
}


# Function to generate sub-network from real network, labeled as 8
graphon8 <- function(u,size){
  #read real data
  RDS_path <- system.file("rds", "email_Eu_core_adj.RDS", package = "localboot")
  P <- readRDS(RDS_path)
  node_samples <- sample(1:NROW(P),size,replace = TRUE)
  P[node_samples,node_samples]
}


#' Generate a Graphon Probability Matrix
#'
#' This function generates a graphon probability matrix based on a specified graphon type. 
#' Users can control the generation process through various parameters.
#'
#' @param size An integer specifying the size of the network.
#' @param graph_num An integer (default is 1) indicating the graphon type to use. 
#'        Acceptable values are from 1 to 6.
#' @param sampling_on_u A logical value determining if uniform sampling should be used for 'u'. 
#'        Defaults to TRUE. If FALSE, a regular sequence from 0 to 1 is used.
#' @param u_input An optional numeric vector that provides specific values for 'u'.
#'        If NULL (default), 'u' is generated based on 'sampling_on_u'.
#'
#' @return A matrix of probabilities is returned. 
#'
#' @examples
#' # Generate a graphon probability matrix of size 100 using graphon setting 1
#' P = generate_graphon(100, 1)
#'
#' @export
generate_graphon <- function(size, graph_num = 1,sampling_on_u=TRUE,u_input=NULL){
  
  #handle user provided u
  if(!is.null(u_input)){u = u_input}else{
    if(sampling_on_u){
      u = stats::runif(n = size)
    }else{
      u = seq(from = 0, to = 1,length.out=size)
    }
  }
  u = sort(u)
  # Use switch to select the graphon
  p_matrix <- switch(as.character(graph_num),
                           "1" = graphon1(u, size),
                           "2" = graphon2(u, size),
                           "3" = graphon3(u, size),
                           "4" = graphon4(u, size),
                           "5" = graphon5(u, size),
                           "6" = graphon6(u, size),
                           "7" = graphon7(u, size),
                           "8" = graphon8(u, size),
                           stop("Invalid graph_num: should be between 1 and 6"))
  return(p_matrix)

}

#' Generate Network Adjacency Matrix from Probability Matrix
#'
#' This function generates a network adjacency matrix from a given probability matrix. 
#' It checks if the input is a valid probability matrix and can produce either a single 
#' network or multiple replicates.
#'
#' @param P A square matrix representing the probability matrix, where each element 
#'     is a probability (between 0 and 1) of an edge between nodes.
#' @param replicate An integer indicating the number of network replicates to generate. 
#'     Defaults to 1.
#' @param symmetric.out A logical value indicating whether the output matrix should be 
#'     symmetric. Defaults to TRUE.
#' @param noloop A logical value indicating whether to include self-loops in the network. 
#'     Defaults to FALSE.
#'
#' @return If `replicate` is 1, returns a single adjacency matrix. If `replicate` is 
#'     greater than 1, returns a list of adjacency matrices. Each matrix is a square 
#'     binary matrix, where 1 indicates the presence of an edge and 0 indicates its absence.
#'
#' @examples
#' P = generate_graphon(100, 1)
#' network = generate_network_P(P, replicate = 1, symmetric.out = TRUE)
#'
#' @export
generate_network_P = function(P, replicate = 1, symmetric.out=TRUE, noloop = FALSE){
  ## Check P
  cond1 = ((all(P>=0))&&(all(P<=1)))
  cond2 = (nrow(P)==ncol(P))
  if (!(cond1&&cond2)){
    stop("* gmodel.P : P is not a valid probability matrix.")
  }
  
  ## Parameter
  n = nrow(P)
  
  ## replicate 1 case
  if (replicate==1){
    tmpmat = matrix(stats::runif(n^2),nrow=n)
    if (symmetric.out){
      tmpmat[lower.tri(tmpmat)] <- t(tmpmat)[lower.tri(t(tmpmat))]
    }
    G = (tmpmat<P)*1
  } else {
    G = list()
    for (i in 1:replicate){
      tmpmat = matrix(stats::runif(n^2),nrow=n)
      if (symmetric.out){
        tmpmat[lower.tri(tmpmat)] <- t(tmpmat)[lower.tri(t(tmpmat))]
      }
      tmpG = 1*(tmpmat<P)
      if (noloop){
        diag(tmpG) = 0
      }
      G[[i]] = tmpG
    }
  }
  ## return output
  return(G)
}
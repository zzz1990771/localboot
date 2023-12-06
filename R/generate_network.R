#this function generate network adj matrix from probability matrix p
generate_network_P = function(P, replicate = 1, symmetric.out=TRUE){
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
    tmpmat = matrix(runif(n^2),nrow=n)
    if (symmetric.out){
      tmpmat[lower.tri(tmpmat)] <- t(tmpmat)[lower.tri(t(tmpmat))]
    }
    G = (tmpmat<P)*1
  } else {
    G = list()
    for (i in 1:replicate){
      tmpmat = matrix(runif(n^2),nrow=n)
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




#generate p matrix from preset graphon function
graphon3 <- function(size){
  
  #param
  return_u = FALSE
  sampling_on_u=TRUE
  u_input=NULL
  
  #handle user provided u
  if(!is.null(u_input)){u = u_input}else{
    if(sampling_on_u){
      u = runif(n = size)
    }else{
      u = seq(from = 0, to = 1,length.out=size)
    }
  }
  u = sort(u)
  
  #graphon function
  smoothness = 8
  p_matrix = matrix(0,nrow=size,ncol=size)
  for(i in 1:size){
    for(j in 1:size){
      p <- cos(pi*(u[i]-u[j]))/smoothness + 0.5
      p_matrix[i,j] <- p    
    }
  }
  
  if(return_u==TRUE){
    return(list(p_matrix=p_matrix,u=u))
  }else{
    return(p_matrix)
  }
  
}
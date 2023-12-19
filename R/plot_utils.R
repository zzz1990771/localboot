#' Plot Adjacency Matrix
#'
#' This function creates a plot of an adjacency matrix, where the matrix is displayed as an image.
#'
#' @param X An adjacency matrix to be plotted.
#' @param ... Additional graphical parameters to pass to `image`.
#'
#' @return Generates a plot.
#'
#' @examples
#' adj_matrix <- matrix(rbinom(100, 1, 0.5), 10, 10)
#' plot_adj(adj_matrix)
#'
#' @export
plot_adj <- function(X,...){
  graphics::image(t(apply(X, 2, rev)),breaks=seq(0, 1, length.out=13),...)
}

#' Plot Probability Matrix
#'
#' This function creates a filled contour plot of a probability matrix using the `viridis` color palette.
#' The plot is created using the `ggplot2` and `viridis` libraries.
#'
#' @param P A probability matrix to be plotted.
#'
#' @return Generates a filled contour plot.
#'
#' @examples
#' P_matrix <- matrix(runif(100), 10, 10)
#' plot_P(P_matrix)
#'
#' @importFrom viridis viridis
#' @importFrom graphics filled.contour
#' @export
plot_P <- function(P){
  #util function  
  rotate <- function(x) t(apply(x, 2, rev))
  graphics::filled.contour(rotate(P),color.palette=viridis::viridis,axes=FALSE,zlim=c(0,1),nlevels=100)
}



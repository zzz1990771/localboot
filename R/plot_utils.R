plot_adj <- function(X,...){
  #hmcols<-colorRampPalette(c("red","white","blue"))(256)
  image(t(apply(X, 2, rev)),breaks=seq(0, 1, length.out=13),...)
}


##plot functions
rotate <- function(x) t(apply(x, 2, rev))
plot_P <- function(P){
  library(viridis)
  library(ggplot2)
  filled.contour(rotate(P),color.palette=viridis,axes=FALSE,zlim=c(0,1),nlevels=100)
}

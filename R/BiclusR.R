#' @title Generating a single bicluster
#' @description Function for generating a single bicluster.
#' @param type a character string indicating the type of bicluster and must be one of "constant", "add", "multi" and "linear" (character).
#' @param a the number of rows in a bicluster (integer).
#' @param b the number of columns in a bicluster (integer).
#' @param M the amplified multiple of value (numeric).
#' @param Mc the amplified multiple of difference for add model (numeric).
#' @param Me the amplified multiple of error (numeric).
#' @details This function can generate four type bicluster. If \code{type} is "\code{constant}", it will output a bicluster with same value. If \code{type} is "\code{add}", it will output a bicluster with a constant difference between each sample. If \code{type} is "\code{multi}", it will output a bicluster with a constant multiple between each sample. If \code{type} is "\code{linear}", it will output a bicluster with a linear relationship between samples.
#' @return a bicluster matrix with dimension (a,b).
#' @examples
#' \dontrun{
#' bictype(type="constant", a = 20, b = 50, M = 5, Me = NA)
#' }
#' @import stats
#' @export
bictype <- function(type, a, b, M, Mc=NA, Me=NA){
  if(type=="constant") x <- matrix(rep(rnorm(1)*M,a*b),nrow=a)
  if(type=="add"){
    x <- matrix(rep(rnorm(b)*M,a),nrow=a,byrow=T)
    for(i in 2:a) x[i,] <- x[i,]+rep(rnorm(1),b)*Mc
  }
  if(type=="multi"){
    x <- matrix(rep(rnorm(b)*M,a),nrow=a,byrow=T)
    for(i in 2:a) x[i,] <- x[i,]*runif(1,0.5,1.5)
  }
  if(type=="linear"){
    x <- matrix(rep(rnorm(b)*M,a),nrow=a,byrow=T)
    for(i in 2:a) x[i,] <- x[i,]*runif(1,0.5,1.5)+rep(rnorm(1),b)*Mc
  }
  
  if(!is.na(Me)) x <- x+matrix(rnorm(a*b)*Me,nrow=a)
  x
}

#' @title Generating a bicluster sample in different structures
#' @description Function for generating a bicluster sample in different structures.
#' @param struct a character string indicating the type of bicluster structure and must be one of "single", "mono.rc" and "chessboard" (character).
#' @param type a character vector indicating the type of each bicluster (vector).
#' @param xdim the dimension of the sample (vector).
#' @param bicdim a list contains two vector naming "row" and "col" which indicate the dimension of each bicluster (list).
#' @param M the amplified multiple of value for each bicluster (vector).
#' @param Mc the amplified multiple of difference for add model for each bicluster (vector).
#' @param Me the amplified multiple of error for each bicluster (vector).
#' @return a bicluster sample in the specified structure.
#' @examples
#' \dontrun{
#' library(pheatmap)
#' type <- rep("constant",3)
#' xdim <- c(200,600)
#' bicdim <- list(row=c(30,60,80),col=c(50,80,120))
#' x <- bicstr("mono.rc",type,xdim,bicdim,M=rep(5,3),rep(2,3))
#' pheatmap(x,show_colnames =F,show_rownames = F, cluster_rows = F, cluster_cols = F)
#' }
#' @import stats
#' @export
bicstr <- function(struct, type, xdim, bicdim,
                   M, Mc=rep(NA,length(type)), Me=rep(NA,length(type))){
  x <- matrix(rnorm(prod(xdim)), nrow=xdim[1])
  
  if(struct=="single") 
    x[1:bicdim$row,1:bicdim$col] <- bictype(type, bicdim$row, bicdim$col, M, Mc, Me)
  
  if(struct=="mono.rc"){
    d <- length(type)
    s <- c(1,1)
    for(i in 1:d){
      x[s[1]:(s[1]+bicdim$row[i]-1),s[2]:(s[2]+bicdim$col[i]-1)] <-
        bictype(type[i], bicdim$row[i], bicdim$col[i], M[i], Mc[i], Me[i])
      s <- s+c(bicdim$row[i],bicdim$col[i])
    }
  }
  
  if(struct=="chessboard"){
    s <- c(1,1)
    t <- 0
    for(i in 1:length(bicdim$row)){
      s[2] <- 1
      for (j in 1:length(bicdim$col)) {
        t <- t+1
        x[s[1]:(s[1]+bicdim$row[i]-1),s[2]:(s[2]+bicdim$col[j]-1)] <-
          bictype(type[t], bicdim$row[i], bicdim$col[j], M[t], Mc[t], Me[t])
        s[2] <- s[2]+bicdim$col[j]
      }
      s[1] <- s[1]+bicdim$row[i]
    }
  }
  x
}
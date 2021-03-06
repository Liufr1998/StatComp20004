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

#' @title evaluate the results from the bicluster method
#' @description Function for evaluating the results from the bicluster method.
#' @param bic.out the biclusters obtained from the bicluster method (list).
#' @param bic.index the true position (row and column) of each bicluster (list).
#' @param bic.type the type of each bicluster (vector).
#' @return a list of information about the matched true group and bicluster.
#' @examples
#' \dontrun{
#' library(biclust)
#' 
#' xdim <- c(300, 600)
#' bicdim <- list(row=c(20,40,70), col=c(50,70,100))
#' bic.index <- list(row=list(1:20,21:60,61:130),
#'                   col=list(1:50,51:120,121:220))
#' type <- rep("constant",3)
#' M <- rep(5,3)
#' Me <- rep(1.5,3)
#' x <- bicstr("mono.rc", type, xdim, bicdim, M, Me=Me)
#' rownames(x) <- paste0("Cell",1:nrow(x))
#' colnames(x) <- paste0("Gene",1:ncol(x))
#' 
#' res <- biclust::biclust(x, method=BCCC())
#' if(res@Number==0) bic.out <- NA
#' if(res@Number!=0){
#'   bic.out <- list()
#'   for(jj in 1:res@Number){
#'     bic <- biclust::bicluster(x, res)[[jj]]
#'     rr <- as.numeric(na.omit(as.numeric(unlist(strsplit(rownames(bic), split="Cell")))))
#'     cc <- as.numeric(na.omit(as.numeric(unlist(strsplit(colnames(bic), split="Gene")))))
#'     bic.out[["row"]][[paste0("clus",jj)]] <- rr
#'     bic.out[["col"]][[paste0("clus",jj)]] <- cc
#'   }
#'  }
#' accuracy(bic.out, bic.index, type)
#' }
#' @import biclust
#' @export
accuracy <- function(bic.out, bic.index, bic.type){
  nclus <- length(bic.out$row)
  ngroup <- length(bic.index$row)
  # same ratio
  clus <- paste0("clus",1:nclus)
  group <- paste0("group",1:ngroup)
  ratio.row <- ratio.col <- ratio <- matrix(nrow=nclus,ncol=ngroup, dimnames=list(clus,group))
  
  for(i in 1:nclus){
    for(j in 1:ngroup){
      ratio.row[i,j] <- length(which(bic.index$row[[j]]%in%bic.out$row[[i]]==T))/length(bic.index$row[[j]])
      ratio.col[i,j] <- length(which(bic.index$col[[j]]%in%bic.out$col[[i]]==T))/length(bic.index$col[[j]])
      ratio[i,j] <- mean(c(ratio.row[i,j], ratio.col[i,j]))
    }
  }
  
  # match group and cluster
  out <- data.frame()
  out[1,c("bictype","id.group","id.clus","ratio.row","ratio.col","average")] <- rep(NA,6)
  for(i in 1:min(ngroup,nclus)){
    ratio2 <- matrix(ratio[clus,group],nrow=length(clus),dimnames=list(clus,group))
    max.ind <- which(ratio2==max(ratio2),arr.ind=T)
    if(nrow(max.ind)>1) max.ind <- max.ind[1,]
    out[i,"id.clus"] <- as.numeric(unlist(strsplit(clus[max.ind[1]], split="clus"))[2])
    out[i,"id.group"] <- as.numeric(unlist(strsplit(group[max.ind[2]], split="group"))[2])
    out[i,"bictype"] <- bic.type[out[i,"id.group"]]
    out[i,"ratio.row"] <- ratio.row[out[i,"id.clus"],out[i,"id.group"]]
    out[i,"ratio.col"] <- ratio.col[out[i,"id.clus"],out[i,"id.group"]]
    out[i,"average"] <- mean(c(out[i,"ratio.row"],out[i,"ratio.col"]))
    
    clus <- clus[-max.ind[1]]
    group <- group[-max.ind[2]]
  }
  out
}
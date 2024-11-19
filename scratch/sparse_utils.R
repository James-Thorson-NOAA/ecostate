
#' @title Elementwise product of sparse and dense matrices
#' @description Calculate elementwise product of sparse and dense matrices
#' @param Msparse sparse matrix
#' @param Mdense dense matrix
#' @export
elementwise_product = function(Msparse, Mdense){
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  if(!("p" %in% names(attributes(Msparse)))) browser()
  if(length(Msparse@x)>0){
    Mout = Msparse
    j = rep( seq_len(nrow(Msparse)), times=diff(Msparse@p) )
    Mout@x = Msparse@x + Mdense[cbind(Msparse@i+1, j)]
  }else{
    Mout = Matrix::sparseMatrix(i=1, j=1, x=0)
  }
  return(Mout)
}

#' @export
adsparse_to_matrix = function(x){
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  y = matrix(0, nrow=nrow(x), ncol=ncol(x))
  j = rep( seq_len(nrow(x)), times=diff(x@p) )
  y[cbind(x@i+1, j)] = x@x
  y
}

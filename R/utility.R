
#' @title Calculate tracers, e.g., trophic level
#'
#' @description Calculate how a tracer propagates through consumption.
#'
#' @inheritParams ecostate
#' @inheritParams compute_nll
#' @inheritParams ecostate_control
#'
#' @param Q_ij Consumption of each prey i by predator j in units biomass.
#' @param tracer_i an indicator matrix specifying the traver value
#'
#' @details
#' Trophic level \eqn{y_i} for each predator \eqn{i} is defined as:
#'
#' \deqn{ \mathbf{y = l Q^* + 1} }
#'
#' where \eqn{\mathbf{Q*}} is the proportion consumption for each predator (column)
#' of different prey (rows).  We identify primary producers as any taxa with no
#' consumption (a column of 0s), and assign them as the first trophic level.
#'
#' More generically, a tracer might be used to track movement of biomass through
#' consumption.  For example, if we have a tracer \eqn{x_i} that is 1 for the 
#' base of the pelagic food chain, and 0 otherwise, then we can calculate 
#' the proportion of pelagic vs. nonpelagic biomass for each taxon:
#'
#' \deqn{ \mathbf{y = l Q^* + x} }
#'
#' This then allows us to separate alternative components of the foodweb.
#'
#' @return
#' The vector \deqn{\mathbf{y_i}} resulting from tracer \deqn{\mathbf{x_i}}
#'
#' @export
compute_tracer <-
function( Q_ij,
          inverse_method = c("Penrose_moore", "Standard"),
          type_i,
          tracer_i = rep(1,nrow(Q_ij)) ){
  # Start
  inverse_method = match.arg(inverse_method)
  assertDouble( tracer_i, lower=0, upper=1, len=nrow(Q_ij), any.missing=FALSE )
  
  # Indicators 
  which_primary = which( type_i=="auto" )
  which_detritus = which( type_i=="detritus" )
  
  # Rescale consumption
  colsums = colSums(Q_ij)
  # diag(vector) doesn't work when length(vector)=1, so using explicit construction
  #B = diag(1/colsums)
  B = matrix(0, nrow=length(colsums), ncol=length(colsums))
  B[cbind(seq_along(colsums),seq_along(colsums))] = 1/colsums
  #inverse_denom = rep(1,nrow(Q_ij)) %*% t(1/colsums)
  Qprime_ij = Q_ij %*% B

  # Identify primary producers
  Qprime_ij[,c(which_primary,which_detritus)] = 0
  
  # Solve for trophic level
  if( inverse_method == "Penrose_moore" ){
    #inverse_IminusQ = corpcor::pseudoinverse( diag(nrow(Qprime_ij)) - Qprime_ij )
    inverse_IminusQ = ginv( diag(nrow(Qprime_ij)) - Qprime_ij )
  }else{
    inverse_IminusQ = solve( diag(nrow(Qprime_ij)) - Qprime_ij )
  }
  x_i = t(tracer_i) %*% inverse_IminusQ
  return( x_i )
}

#' @title Penrose-Moore pseudoinverse
#' @description Extend \code{MASS:ginv} to work with RTMB
#' @param x Matrix used to compute pseudoinverse
#' @importFrom MASS ginv
ginv <- RTMB::ADjoint(function(x) {
    n <- sqrt(length(X))
    dim(X) <- c(n,n)
    MASS::ginv(X)
}, function(X,Y,dY) {
    n <- sqrt(length(X))
    dim(Y) <- dim(dY) <- c(n,n)
    -t(Y)%*%dY%*%t(Y)
}, name="ginv")

#' @title Dirichlet-multinomial
#' @description Allows data-weighting as parameter
#' @param x numeric vector of observations across categories
#' @param prob numeric vector of category probabilities
#' @param ln_theta logit-ratio of effective and input sample size
#' @param log whether to return the log-probability or not
#' @examples
#' library(RTMB)
#' prob = rep(0.1,10)
#' x = rmultinom( n=1, prob=prob, size=20 )[,1]
#' f = function( ln_theta ) ddirmult(x, prob, ln_theta)
#' f( 0 )
#' F = MakeTape(f, 0)
#' F$jacfun()(0)
#'
#' @return
#' The log-likelihood resulting from the Dirichlet-multinomial distribution
#'
#'
#' @export
ddirmult <-
function( x,
          prob,
          ln_theta,
          log=TRUE ){

  # Pre-processing
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  Ntotal = sum(x)
  p_exp = prob / sum(prob)
  p_obs = x / Ntotal
  dirichlet_Parm = exp(ln_theta) * Ntotal
  logres = 0.0

  # https://github.com/nmfs-stock-synthesis/stock-synthesis/blob/main/SS_objfunc.tpl#L306-L314
  # https://github.com/James-Thorson/CCSRA/blob/master/inst/executables/CCSRA_v8.cpp#L237-L242
  # https://www.sciencedirect.com/science/article/pii/S0165783620303696

  # 1st term -- integration constant that could be dropped
  logres = logres + lgamma( Ntotal+1 )
  for( index in seq_along(x) ){
    logres = logres - lgamma( Ntotal*p_obs[index] + 1.0 )
  }

  # 2nd term in formula
  logres = logres + lgamma( dirichlet_Parm ) - lgamma( Ntotal+dirichlet_Parm )

  # Summation in 3rd term
  for( index in seq_along(x) ){
    logres = logres + lgamma( Ntotal*p_obs[index] + dirichlet_Parm*p_exp[index] )
    logres = logres - lgamma( dirichlet_Parm * p_exp[index] )
  }
  if(log){ return(logres) }else{ return(exp(logres)) }
}


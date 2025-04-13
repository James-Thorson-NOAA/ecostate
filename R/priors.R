
#' EcoState priors
#'
#' @param priors Either: 
#'  \itemize{
#'    \item A list of sampling statements (e.g., `logtau_i[] ~ dnorm(mean = log(0.5), sd = 1)`) for named parameters. 
#'    See "Details" for a list of parameters and their descriptions.
#'    \item A function that takes as input the list of
#'        parameters \code{out$obj$env$parList()} where \code{out} is the output from
#'        \code{ecostate()}, and returns a numeric vector
#'        where the sum is the log-prior probability.  For example
#'        \code{log_prior = function(p) dnorm( p$logq_i[1], mean=0, sd=0.1, log=TRUE)}
#'        specifies a lognormal prior probability for the catchability coefficient
#'        for the first \code{taxa} with logmean of zero and logsd of 0.1
#'  }
#' @param p List of parameters
#'
#' @return Log-prior density (numeric)
#' @export
#'
#' @details
#' Priors can be specified as a list of sampling statements, with parameter names on the left-hand side and density 
#' functions and their named arguments on the right-hand side. SEM parameters should be referred to by their names,
#' while all other parameters are indexed (e.g., by taxa, year, or stanza group) and should be specified with brackets.
#' E.g., the following list specifies priors on an SEM parameter named \code{rho_X}, the biomass process errors for all time
#' steps for taxa Y \code{epsilon_ti[,"Y"]}, and the ecotrophic efficiency for all taxa:
#' \preformatted{
#'  priors <- list(
#'    rho_X ~ dnorm(mean = 0, sd = 0.5), 
#'    epsilon_ti[,"Y"] ~ dnorm(mean = 0, sd = 0.1), 
#'    EE_i[] ~ dnorm(mean = 1, sd = 0.1)
#'  )
#' }
#' 
evaluate_prior <- function(priors, p) {
  
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  
  if (length(p$beta) > 0) { 
    for (i in seq_along(p$beta)) {
      if (isTRUE(names(p$beta[i]) %in% names(p))) stop(paste("Parameter name", names(p$beta[i]), "is reserved"))
      p[[names(p$beta[i])]] <- p$beta[i]
    }
    p <- p[names(p) != "beta"]
  }
  
  if (inherits(priors, "function")) {
    
    logp <- priors(p)
    
  } else if (inherits(priors, "list")) {
    
    logp <- 0
    
    for (i in seq_along(priors)) {
      
      # Parse prior density function
      rhs <- priors[[i]][[3]]
      dens <- as.character(rhs[[1]])
      
      # Unpack density arguments
      args <- as.list(formals(dens))
      args <- modifyList(args, as.list(rhs)[-1])
      args$x <- priors[[i]][[2]]
      args$log <- TRUE
      
      # Evaluate prior density
      logp_i <- try(with(p, sum(do.call(dens, args))), silent = TRUE)
      if (inherits(logp_i, "try-error") | is.na(logp_i)) stop(paste("Problem with prior", deparse(args$x)))
      logp <- logp + logp_i
    }
    
  }
  
  return(logp)
  
}

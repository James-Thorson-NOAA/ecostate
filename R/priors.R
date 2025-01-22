
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
#' @param taxa Character vector of taxa included in mode
#' @param years Integer-vector of years included in model
#' @param stanza_groups Character vector of unique multi-stanza groups
#' @param sem SEM model data frame
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
evaluate_prior <- function(priors, p, taxa, years, stanza_groups, sem = "") {
  
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  
  if (class(priors) == "function") {
    
    logp <- priors(p)
    
  } else if (class(priors) == "list") {
    
    logp <- 0
    
    for (i in seq_along(priors)) {
      
      # Parse formula
      lhs <- priors[[i]][[2]]
      rhs <- priors[[i]][[3]]
      dens <- as.character(rhs[[1]])
      
      # Unpack density arguments
      args <- as.list(formals(dens))
      args <- modifyList(args, as.list(rhs)[-1])
      args$log <- TRUE
      
      # Match LHS to SEM parameter names
      if (class(sem) == "data.frame") {
        
        if (as.character(as.list(lhs)[[1]]) %in% sem$name) {
          
          args$x <- p$beta[unique(sem$parameter[which(sem$name == as.character(lhs))])]
          logp <- logp + do.call(dens, args)
          
        }
        
      }
      
      # Match indexed parameters
      if (as.character(as.list(lhs)[[1]]) == "[") {
        
        lhs <- as.list(lhs)
        par_i <- as.character(lhs[[2]])
        
        if (!(par_i %in% names(p))) stop(paste0("Parameter name `", par_i, "` not recognized"))  
        
        if (grepl("_i", par_i)) {
          
          taxa_i <- if (all(as.character(lhs[[3]]) == "")) seq_along(taxa) else do.call("match", list(x = lhs[[3]], table = taxa))
          args$x <- p[[par_i]][taxa_i]
          logp <- logp + sum(do.call(dens, args))
          
        } else if (grepl("_ti", par_i)) {
          
          year_i <- if (all(as.character(lhs[[3]]) == "")) seq_along(years) else do.call("match", list(x = lhs[[3]], table = years))
          taxa_i <- if (all(as.character(lhs[[4]]) == "")) seq_along(taxa) else do.call("match", list(x = lhs[[4]], table = taxa))
          args$x <- p[[par_i]][year_i, taxa_i]
          logp <- logp + sum(do.call(dens, args))
          
        } else if (grepl("_g2", par_i)) {
          
          sg_i <- if (all(as.character(lhs[[3]]) == "")) seq_along(stanza_groups) else do.call("match", list(x = lhs[[3]], table = stanza_groups))
          args$x <- p[[par_i]][sg_i]
          logp <- logp + sum(do.call(dens, args))
          
        } else if (grepl("_tg2", par_i)) {
          
          year_i <- if (all(as.character(lhs[[3]]) == "")) seq_along(years) else do.call("match", list(x = lhs[[3]], table = years))
          sg_i <- if (all(as.character(lhs[[4]]) == "")) seq_along(stanza_groups) else do.call("match", list(x = lhs[[4]], table = stanza_groups))
          args$x <- p[[par_i]][year_i, sg_i]
          logp <- logp + sum(do.call(dens, args))
          
        } else if (grepl("_z", par_i)) {
          
          # To implement
          
        }
        
      }
      
    }
    
  }
  
  if (is.na(logp)) stop("Problem with specified prior(s)")
  
  return(logp)
  
}
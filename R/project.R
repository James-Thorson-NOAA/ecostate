
#' @title Conditional simulation from a GMRF
#'
#' @description
#' Generates samples from a Gaussian Markov random field (GMRF) conditional upon
#' fixed values for some elements. Modified from package \code{tinyVAST}.
#'
#' @param Q precision for a zero-centered GMRF.
#' @param observed_idx integer vector listing rows of \code{Q} corresponding to
#'        fixed measurements
#' @param x_obs numeric vector with fixed values for indices \code{observed_idx}
#' @param nsim integer listing number of simulated values
#' @param what Whether to simulate from the conditional GMRF, or predict the mean
#'        and precision
#'
#' @return
#' A matrix with \code{nsim} columns and a row for every row of \code{Q} not in
#' \code{observed_idx}, with simulations for those rows
#'
#' @export
conditional_gmrf <-
  function( Q,
            observed_idx,
            x_obs,
            nsim = 1,
            what = c("simulate","predict") ){
    
    # Required libraries
    #library(Matrix)
    what = match.arg(what)
    
    # Error checks
    if( !all(observed_idx %in% seq_len(nrow(Q))) ){
      stop("Check `observed_idx` in `conditional_gmrf")
    }
    if( length(observed_idx) != length(x_obs) ){
      stop("Check length of `observed_idx` and `x_obs`")
    }
    if( any(is.na(x_obs)) ){
      stop("`x_obs` cannot include NA values")
    }
    
    # Calculate conditional mean and variance
    predict_conditional_gmrf <- function(Q, observed_idx, x_obs) {
      all_idx <- seq_len(nrow(Q))
      unobserved_idx <- setdiff(all_idx, observed_idx)
      
      # Partition Q
      #Q_oo <- Q[observed_idx, observed_idx, drop = FALSE]
      Q_ou <- Q[observed_idx, unobserved_idx, drop = FALSE]
      Q_uo <- Matrix::t(Q_ou)
      Q_uu <- Q[unobserved_idx, unobserved_idx, drop = FALSE]
      
      # Compute conditional mean and covariance
      #mu_cond <- -Q_uu_inv %*% Q_uo %*% x_obs
      mu_cond <- -1 * Matrix::solve(Q_uu, Q_uo %*% x_obs)
      
      out = list( mean = as.vector(mu_cond),
                  Q_uu = Q_uu,
                  unobserved_idx = unobserved_idx )
      return(out)
    }
    
    res <- predict_conditional_gmrf(Q, observed_idx, x_obs )
    if( what == "predict" ){
      return(res)
    }else{
      y = c(res$mean + rgmrf0(n = nsim, Q = res$Q_uu))
      return(y)
    }
  }

# Unexported
# Update first n values in atomics in list x with those in list y
substitute_pars <- function(x, y) {
  
  ind <- intersect(names(x), names(y))
  
  for (i in ind) {
    
    yd <- dim(y[[i]])
    
    if (is.null(yd)) {
      # Vectors
      x[[i]][seq_len(length(y[[i]]))] <- y[[i]]
    } else {
      # Matrices / arrays
      yind <- lapply(yd, seq_len)
      args <- c(list(x[[i]]), yind, list(value = y[[i]]))
      x[[i]] <- do.call("[<-", args)
    }
  }
  return(x)
}

#' @title Project ecostate model to future times (EXPERIMENTAL)
#'
#' @description
#' Projects a fitted model forward in time, optionally 
#' conditioning on covariates and/or catches.
#'
#' @param model fitted model from \code{ecostate}
#' @param nsim Number simulations
#' @param extra_years a vector of extra years, matching values in \code{Frate} and 
#'        \code{covariates}
#' @param Frate matrix with rownames corresponding to \code{extra_years} and colnames
#'        corresponding to \code{model$internal$taxa} giving fishing mortality in future
#'        years. Any missing years or taxa will be assigned an F rate of 0.
#' @param covariates matrix with rownames corresponding to \code{extra_years}
#'        and columns corresponding to SEM covariates. Covariates corresponding to 
#'        missing covariate columns and/or NA values within included covariate columns
#'        will be estimated conditional on included covariates.
#' @param future_var logical indicating whether to simulate future process errors
#'        from GMRF, or just compute the predictive mean. If \code{covariates} are
#'        provided for future years, process errors will be simulated conditional
#'        on them.
#' @param past_var One of "fixed", "conditional", or "unconditional", where "fixed" holds
#'        past process errors at their estimates, "conditional" re-simulates them
#'        conditional on past covariates in the model, and "unconditional" re-simulates
#'        both past process errors and covariates from their joint distribution. 
#'        "conditional" and "unconditional" past process errors thus change the boundary 
#'        conditions of the projection.
#' @param parm_var logical indicating whether to re-sample fixed effects from their
#'        predictive distribution, thus changing the GMRF for future process errors
#' @param obj_new Optional; a re-built model object with parameter blocks covering
#'        the full range of years in the simulation. Intended for internal use.
#'
#' @return
#' A list (or list of lists, if nsim > 1) of estimated parameters and 
#' derived quantities with structure resembling \code{model$internal$parhat}
#' 
#' @export
project <- function(
    model, 
    nsim = 1,
    extra_years,
    Frate,
    covariates, 
    future_var = TRUE, 
    past_var = c("fixed", "conditional", "unconditional"), 
    parm_var = FALSE, 
    obj_new
) {
  
  past_var <- match.arg(past_var)
  
  taxa <- model$internal$taxa
  stgroups <- model$internal$settings$unique_stanza_groups
  parlist <- model$internal$parhat
  years_all <- c(model$internal$years, extra_years)
  
  if (!identical(min(years_all):max(years_all), years_all)) {
    stop("extra_years must be an integer sequence starting from the year following the last year in the model.")
  }
  
  if (!is.null(model$internal$extra_years)) {
    stop("`project` is mutually exclusive with `future` argument")
  }

  # For unconditional past variance, resample all process errors and covariates
  if (isFALSE(parm_var) & past_var == "unconditional") {
    parsim <- model$simulator(parlist, simulate_random = TRUE)
    parlist <- substitute_pars(parlist, parsim)
  }
  
  # For conditional, drop process errors from parlist 
  # to be imputed later based on covariates
  if (past_var == "conditional") {
    parlist$epsilon_ti[!is.na(model$tmb_inputs$map$epsilon_ti)] <- NA
    parlist$nu_ti[!is.na(model$tmb_inputs$map$nu_ti)] <- NA
    parlist$nu_tij[!is.na(model$tmb_inputs$map$nu_tij)] <- NA
    parlist$phi_tg2[!is.na(model$tmb_inputs$map$phi_tg2)] <- NA
  }
  
  if (isTRUE(parm_var)) message("resmapling of fixed effects not yet implemented")
  
  # Populate / format logF_ti matrix
  if (!missing(Frate)) {
    
    logF_new <- matrix(
      -Inf, nrow = length(extra_years), ncol = length(taxa), 
      dimnames = list(extra_years, taxa)
    )
    
    if (!all(colnames(Frate) %in% taxa)) {
      message(paste(
        "ignoring Frate for taxa not included in the model:",
        paste(colnames(Frate)[!(colnames(Frate) %in% taxa)], collapse = ", ")
      ))
      Frate <- Frate[,colnames(Frate) %in% taxa, drop = FALSE]
    }
    
    if (!all(rownames(Frate) %in% as.character(extra_years))) {
      message("ignoring years in Frate matrix not in extra_years")
      Frate <- Frate[rownames(Frate) %in% as.character(extra_years),,drop = FALSE]
    }
    
    for (i in seq_along(colnames(Frate))) {
      logF_new[rownames(Frate),colnames(Frate)[i]] <- log(Frate[,colnames(Frate)[i]])
    }
  
    parlist$logF_ti <- rbind(parlist$logF_ti, logF_new)
    rownames(parlist$logF_ti) <- years_all
      
  }
  
  # Populate / format covariates matrix
  cov_hind <- parlist$covariates
  
  cov_new <- matrix(
    NA, nrow = length(extra_years), ncol = ncol(cov_hind), 
    dimnames = list(extra_years, colnames(cov_hind))
  )
  
  if (!missing(covariates)) {
    
    if (!all(colnames(covariates) %in% colnames(cov_hind))) {
      message(paste(
        "ignoring covariates not included in the model:",
        paste(colnames(covariates)[!(colnames(covariates) %in% colnames(cov_hind))], collapse = ", ")
      ))
      covariates <- covariates[,colnames(covariates) %in% colnames(cov_hind), drop = FALSE]
    }
    
    if (!all(rownames(covariates) %in% as.character(extra_years))) {
      message("ignoring years in covariates matrix not in extra_years")
      covariates <- covariates[rownames(covariates) %in% as.character(extra_years),,drop = FALSE]
    }
    
    for (i in seq_along(colnames(covariates))) {
      cov_new[rownames(covariates),colnames(covariates)[i]] <- covariates[,colnames(covariates)[i]]
    }
    
  } 
  
  cov_full <- rbind(cov_hind, cov_new)
  
  # Rebuild object, if not provided
  if (missing(obj_new)) {
    # Control settings 
    ctrl_new <- model$internal$control
    ctrl_new$nlminb_loops <- ctrl_new$newton_loops <- 0
    ctrl_new$getsd <- FALSE
    ctrl_new$map <- ctrl_new$tmb_par <- NULL
    
    # Some objects needed to update the model are not returned by
    # ecostate(), or may not always evaluate correctly (e.g. values from call()). 
    # Use update() to recycle most arguments by calling in the parent env., 
    # but pull specific values to use from the function environment
    new_call <- bquote(update(
      object = .(model), 
      years = .(years_all), 
      covariates = .(cov_full), 
      control = .(ctrl_new)
    ))
    
    # Build
    obj_new <- eval(new_call, envir = parent.frame())
  }
  
  # Replace values in years prior to projection with fitted estimates
  par_new <- obj_new$internal$parhat 
  par_new <- substitute_pars(par_new, parlist)
  
  # Predict random effects in future years conditional on fitted estimates
  yr_ind <- seq_along(model$internal$years)
  
  if (!is.null(model$sem)) {
    
    sem_path <- eval(model$call$sem, envir = parent.frame())
    
    # Store back to model call, to avoid problems with recursion when nsim > 1
    model$call$sem <- bquote(.(sem_path))
    
    sem <- parse_ecostate_sem(
      sem_path, covariates = cov_full, taxa = taxa, 
      years = years_all, fit_eps = c(), fit_nu = c(), 
      settings = model$internal$settings,
      control = model$internal$control
    )$model
    
    # SEM precision matrix
    variables <- unique(c(sem$first, sem$second))
    sem_mat <- make_matrices(par_new$beta, sem, years_all, variables)
    Q <- t(sem_mat$IminusP_kk) %*% sem_mat$invV_kk %*% sem_mat$IminusP_kk
    
    # Matrix to bundle process errors and covariates
    Xit <- matrix(NA, nrow = length(years_all), ncol = length(variables))
    colnames(Xit) <- variables
    
    # Subtract estimated covariate means
    if (!is.null(dim(cov_full))) {
      Xit[,colnames(cov_full)] <- cov_full
      Xit[,colnames(cov_full)] <- sweep(Xit[,colnames(cov_full), drop = FALSE], 2, par_new$mu) 
    }
        
    # Pull out process error values from matrices
    for (i in seq_len(ncol(Xit))) {
      
      if (gsub("eps_", "", colnames(Xit)[i]) %in% taxa) {
        Xit[yr_ind,i] <- parlist$epsilon_ti[,which(taxa %in% gsub("eps_", "", colnames(Xit)[i]))]
      } else if (grepl("nu_", colnames(Xit)[i])) {
        if (gsub("nu_", "", colnames(Xit)[i]) %in% taxa) {
          Xit[yr_ind,i] <- parlist$nu_ti[,which(taxa %in% gsub("nu_", "", colnames(Xit)[i]))]
        } else if (all(strsplit(gsub("nu_", "", colnames(Xit)[i]), ":")[[1]] %in% taxa)) {
          pred_prey <- strsplit(gsub("nu_", "", colnames(Xit)[i]), ":")[[1]]
          Xit[yr_ind,i] <- parlist$nu_tij[, pred_prey[1], pred_prey[2]]
        } 
      } else if (gsub("phi_", "", colnames(Xit)[i]) %in% stgroups) {
        Xit[yr_ind,i] <- parlist$phi_tg2[,which(stgroups %in% gsub("phi_", "", colnames(Xit)[i]))]
      }
    }
    
    # Stack columnwise
    Xvec <- c(Xit)
    
    # Treat non-NA indices as fixed and condition the GMRF on them
    X_fixed <- which(!is.na(Xvec))
    
    # Conditional simulation from GMRF
    Xit[is.na(Xit)] <- conditional_gmrf(
      Q, 
      observed_idx = X_fixed, 
      x_obs = Xvec[X_fixed], nsim = 1, 
      what = ifelse(isTRUE(future_var), "simulate", "predict")
    )
    
    # Add back in estimated covariate means
    if (!is.null(dim(cov_full))) {
      Xit[,colnames(cov_full)] <- sweep(Xit[,colnames(cov_full), drop = FALSE], 2, par_new$mu, FUN = "+")
    }
    
    # Unpack samples from GMRF
    for (i in seq_len(ncol(Xit))) {
      if (gsub("eps_", "", colnames(Xit)[i]) %in% taxa) {
        par_new$epsilon_ti[, which(taxa %in% gsub("eps_", "", colnames(Xit)[i]))] <- Xit[,i]
      } else if (grepl("nu_", colnames(Xit)[i])) {
        if (gsub("nu_", "", colnames(Xit)[i]) %in% taxa) {
          par_new$nu_ti[, which(taxa %in% gsub("nu_", "", colnames(Xit)[i]))] <- Xit[,i]
        } else if (all(strsplit(gsub("nu_", "", colnames(Xit)[i]), ":")[[1]] %in% taxa)) {
          pred_prey <- strsplit(gsub("nu_", "", colnames(Xit)[i]), ":")[[1]]
          par_new$nu_tij[, pred_prey[1], pred_prey[2]] <- Xit[,i]
        } 
      } else if (gsub("phi_", "", colnames(Xit)[i]) %in% stgroups) {
        par_new$phi_tg2[, which(stgroups %in% gsub("phi_", "", colnames(Xit)[i]))] <- Xit[,i]
      } else if (colnames(Xit)[i] %in% colnames(par_new$covariates)) {
        par_new$covariates[,which(colnames(par_new$covariates) == colnames(Xit)[i])] <- Xit[,i]
      }
    }
    
  } else if (isTRUE(future_var)) {
    
    # If model does not use SEM,
    # implement i.i.d normal simulation for process errors
    
    for (i in seq_along(taxa)) {
      
      if (!is.na(par_new$logtau_i[taxa[i]])) {
        par_new$epsilon_ti[,taxa[i]] <- rnorm(length(years_all), mean = 0, sd = exp(par_new$logtau_i[taxa[i]]))
      }
      
      if (!is.na(par_new$logsigma_i[taxa[i]])) {
        par_new$nu_ti[,taxa[i]] <- rnorm(length(years_all), mean = 0, sd = exp(par_new$logsigma_i[taxa[i]]))
      }
      
    }
    
    for (i in seq_along(stgroups)) {
      
      if (!is.na(par_new$logpsi_g2[stgroups[i]])) {
        par_new$phi_tg2[,stgroups[i]] <- rnorm(length(years_all), mean = 0, sd = exp(par_new$logpsi_g2[stgroups[i]]))
      }
      
    }
    
  }
  
  # Simulate
  # Random is false because random simulation handled above
  if (nsim == 1) {
    sim <- obj_new$simulator(par_new, simulate_random = FALSE)
  } else {
    # Call project() recursively but provide obj_new to avoid rebuilding
    sim <- vector("list", length = nsim)
    for (i in seq_len(nsim)) {
      sim[[i]] <- project(
        model = model, obj_new = obj_new, nsim = 1,
        extra_years = extra_years, Frate = Frate, covariates = covariates, 
        future_var = future_var, past_var = past_var, parm_var = parm_var
      )
    }
  }
  
  return(sim)
  
}


parse_ecostate_sem <- function(sem, covariates, taxa, years, fit_eps, fit_nu, settings, control) {

  read_model <- utils::getFromNamespace("read_model", "dsem")
  
  if(isTRUE((length(fit_eps) > 0 | length(fit_nu > 0) | length(settings$fit_phi) > 0))) {
    warning("arguments fit_eps, fit_nu, and settings$fit_phi are ignored if using SEM")
  }
  
  if(isTRUE(control$process_error != "epsilon")) {
    warning("SEM forces innovation ('epsilon') process errors")
    control$process_error <<- "epsilon"
  }
  
  paths = paste(scan(
    text = sem, 
    what = list(path = "", lag = 1, par = "", start = 1, dump = ""), 
    sep = ",", strip.white = TRUE, 
    comment.char = "#", fill = TRUE, quiet = TRUE
  )$path, collapse = " ")
  
  sem_vars <- proc_vars <- unique(regmatches(paths, gregexpr("(eps_|nu_|phi_)[^\\s,]+", paths, perl = TRUE))[[1]])
  
  if (!is.null(covariates)) {
    
    if (isFALSE(nrow(covariates) == length(years))) {
      stop("matrix of SEM covariates must have a row for each year")
    }
    
    if (isFALSE(ncol(covariates) == length(unique(colnames(covariates))))) {
      stop("Please check `colnames(covariates)` to confirm that all variables (columns) have a unique name")
    }
    
    sem_vars <- c(colnames(covariates), sem_vars)
    
  }
  
  sem <- read_model(sem, times = years, variables = sem_vars, quiet = TRUE)
  
  if (any(grepl(" ", taxa))) {
    stop("If using SEM, taxa (and stanza) names must not have spaces")
  }
  
  if (length(settings$stanza_groups) > 0) {
    if (any(grepl(" ", settings$unique_stanza_groups))) {
      stop("If using SEM, taxa (and stanza) names must not have spaces")
    }
    for (i in seq_along(proc_vars)) {
      if (gsub("eps_|nu_", "", proc_vars[i]) %in% settings$unique_stanza_groups) {
        if (!(gsub("eps_|nu_", "", proc_vars[i]) %in% taxa)) {
          stop("epsilon and nu process errors must be specified for taxa, not stanza groups")
        }
      }
    }
  }
  
  if (any((sem$direction == 2) & (sem$lag != 0))) {
    stop("All two-headed arrows should have lag=0")
  }
  
  if (!all(c(sem$first, sem$second) %in% sem_vars)) {
    stop("Some variable(s) in `sem` are not in `covariates`, or process errors are not correctly denoted")
  }
  
  # Unique parameters to be estimated
  beta <- unname(c(tapply(sem$start[sem$parameter != 0], sem$parameter[sem$parameter != 0], unique)))
  if (isTRUE(class(beta) == "list")) stop ("Differing starting values provided for SEM parameter(s) with equality constraint")
  beta[is.na(beta)] <- sem$start[is.na(sem$start)] <- 0.1
  beta <- setNames(beta, unique(as.character(na.omit(sem$name))))
  
  return(list(model = sem, beta = beta, proc_vars = proc_vars))
  
}
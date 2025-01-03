
library("tidyverse")
library("checkmate")
library("Matrix")
library("dsem")
library("RTMB")
library("ecostate")

data("eastern_bering_sea", package = "ecostate")
ebs <- eastern_bering_sea

# Set up inputs for ecostate() ----------------------------

# Data
taxa <- names(ebs$P_over_B)
years <- sort(unique(ebs$Survey$Year))
catch <- ebs$Catch
biomass <- ebs$Survey
DC <- ebs$Diet_proportions
PB <- ebs$P_over_B
QB <- ebs$Q_over_B
type <- sapply(taxa, switch, "Detritus" = "detritus", "Chloro" = "auto", "hetero")

# Initial equilibrium biomass and ecotrophic efficiencies
EE <- B <- U <- setNames(rep(NA, length(taxa)), taxa)
B[] <- aggregate(ebs$Survey$Mass, list(Taxon = factor(ebs$Survey$Taxon, levels = taxa)), mean, drop = FALSE)[,2]
EE[] <- 1
U[] <- 0.2 # (?)

# Default vulnerability, except producers (?)
X <- array(2, dim = c(length(taxa), length(taxa)))
dimnames(X) <- list(taxa, taxa)
X[,'Chloro'] <- 91

# Define priors
log_prior <- function(p){
  # Prior on process-error log-SD to stabilize model
  logp = sum(dnorm( p$logtau_i, mean=log(0.2), sd=1, log=TRUE ), na.rm=TRUE)
}

# Taxa to include
taxa <- c( "NFS", "Pollock", "Copepod", "Chloro", "Krill" )

# Parameters to estimate
fit_Q <- c("Pollock", "Copepod", "Chloro", "Krill")
fit_B0 <- c("Pollock", "NFS")
fit_B <- c("Cod", "NFS")  
fit_eps <- "Pollock"

# DSEM paths
sem <- "
  cold_pool -> eps_Pollock, 1, beta_cp_plk, 0,
  eps_Pollock -> eps_Pollock, 1, rho_plk, 0,
  cold_pool <-> cold_pool, 0, sigma_cp, 0.1
  eps_Pollock <-> eps_Pollock, 0, sigma_plk, 0.1
"

# DSEM additional covariates
covariates <- structure(
  c(0.349, 0.278, 0.481, 0.491, 0.462, 0.3, 0.381, 0.181, 
    0.298, 0.379, 0.57, 0.223, 0.379, 0.59, 0.373, 0.56, 0.27, 0.708, 
    0.515, 0.424, 0.379, 0.292, 0.379, 0.351, 0.602, 0.682, 0.69, 
    0.7, 0.702, 0.544, 0.722, 0.615, 0.335, 0.355, 0.27, 0.475, 0.055, 
    0.187, 0.535, 0.295, 0.475, 0.423), 
  dim = c(42L, 1L), 
  dimnames = list(NULL, "cold_pool")
  )

# Other values
agecomp = list()
weight = list()
fit_nu <- vector()
fit_EE <- vector()
fit_PB <- vector()
fit_QB <- vector()
settings <- stanza_settings(taxa=taxa)
control <- ecostate_control(n_steps = 20, start_tau = 0.01, tmbad.sparse_hessian_compress = 0, silent = FALSE)

# Function to make DSEM precision matrix ------------------

make_matrices <-
  function( beta_p,
            model,
            times,
            variables ){
    
    model <- as.data.frame(model)
    
    if(missing(beta_p)){
      model_unique = model[match(unique(model$parameter),model$parameter),]
      beta_p =  as.numeric(model_unique$start)
    }
    
    # Loop through paths
    P_kk = drop0(sparseMatrix( i=1, j=1, x=0, dims=rep(length(variables)*length(times),2) ))   # Make with a zero
    #P_kk = AD(P_kk)
    G_kk = (P_kk)
    for( i in seq_len(nrow(model)) ){
      lag = as.numeric(model[i,2])
      L_tt = sparseMatrix( i = seq(lag+1,length(times)),
                           j = seq(1,length(times)-lag),
                           x = 1,
                           dims = rep(length(times),2) )
      
      P_jj = sparseMatrix( i = match(model[i,'second'],variables),
                           j = match(model[i,'first'],variables),
                           x = 1,
                           dims = rep(length(variables),2) )
      
      tmp_kk = (kronecker(P_jj, L_tt))
      
      # Assemble
      if(abs(as.numeric(model[i,'direction']))==1){
        P_kk = P_kk + beta_p[as.integer(model$parameter[i])] * tmp_kk # AD(tmp_kk)
      }else{
        G_kk = G_kk + beta_p[as.integer(model$parameter[i])] * tmp_kk # AD(tmp_kk)
      }
    }
    
    # Diagonal component
    I_kk = Diagonal(nrow(P_kk))
    
    # Assemble
    IminusP_kk = I_kk - P_kk
    invV_kk = AD(G_kk)
    invV_kk@x = 1 / G_kk@x^2
    Q_kk = t(IminusP_kk) %*% invV_kk %*% IminusP_kk
    
    out = list(
      "P_kk" = P_kk,
      "G_kk" = G_kk,
      "invV_kk" = invV_kk,
      "IminusP_kk" = IminusP_kk,
      "Q_kk" = Q_kk
    )
    return(out)
  }


# Run -----------------------------------------------------

# Current ecostate implementation, with fit_eps
# Running as a check that I don't break compatability with fit_eps syntax...
run_eps <- ecostate::ecostate(
  taxa = taxa, years = years, catch = catch, biomass = biomass, 
  PB = PB, QB = QB, B = B, DC = DC, EE = EE, X = X, 
  type = type, U = U, fit_Q = fit_Q, fit_B0 = fit_B0, fit_B = fit_B, fit_eps = fit_eps, 
  settings = settings, control = control
)

# Expose internal ecostate and RTMB functions
# Shortcut to avoid building ecostate package
library("pacman")
attach(loadNamespace("TMB"), name = "TMB_all")
attach(loadNamespace("RTMB"), name = "RTMB_all")
attach(loadNamespace("ecostate"), name = "ecostate_all")

# Source dev versions
source("R/compute_nll.R")
source("R/ecostate.R")

# Run again with sourced version
run_eps2 <- ecostate(
  taxa = taxa, years = years, catch = catch, biomass = biomass, 
  PB = PB, QB = QB, B = B, DC = DC, EE = EE, X = X, 
  type = type, U = U, fit_Q = fit_Q, fit_B0 = fit_B0, fit_B = fit_B, fit_eps = fit_eps, 
  settings = settings, control = control
)

# Run with SEM
run_sem <- ecostate(
  taxa = taxa, years = years, catch = catch, biomass = biomass, 
  PB = PB, QB = QB, B = B, DC = DC, EE = EE, X = X, 
  type = type, U = U, fit_Q = fit_Q, fit_B0 = fit_B0, fit_B = fit_B, 
  sem = sem, covariates = covariates,
  settings = settings, control = control
)

# Compare parameters --------------------------------------

drop <- c(setdiff(names(run_sem$opt$par), names(run_eps$opt$par)), setdiff(names(run_eps$opt$par), names(run_sem$opt$par)))

par_shared <- cbind(
  eps_main = run_eps$opt$par[!(names(run_eps$opt$par) %in% drop)],
  eps_dev = run_eps$opt$par[!(names(run_eps$opt$par) %in% drop)],
  sem = run_sem$opt$par[!(names(run_sem$opt$par) %in% drop)]
  )

cor(par_shared)

sem <- make_dsem_ram(sem, years, c("cold_pool", "eps_Pollock"))$model
c(sem_coef <- setNames(run_sem$opt$par[4:7], sem[,"name"]))
image(make_matrices(sem_coef, sem, years, c("cold_pool", "eps_Pollock"))$Q_kk)

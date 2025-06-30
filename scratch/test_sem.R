
library("tidyverse")
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
log_prior = list(
  rho_plk ~ dnorm(mean = 0, sd = 0.25),
  sigma_plk ~ dlnorm(meanlog = -2, sdlog = 0.5)
)

# Taxa to include
taxa <- c( "NFS", "Pollock", "Copepod", "Chloro", "Krill", "Arrowtooth", "Cod")

# Parameters to estimate
fit_Q <- c("Pollock", "Copepod", "Chloro", "Krill")
fit_B0 <- c("Pollock", "Arrowtooth", "Cod", "NFS")
fit_EE <- vector()
fit_B <- vector()  #c("Cod", "NFS")

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
fit_PB <- vector()
fit_QB <- vector()
settings <- stanza_settings(taxa=taxa)
control <- ecostate_control(n_steps = 20, tmbad.sparse_hessian_compress = 0)

# Simulation testing, epsilon -----------------------------

# DSEM paths
sem <- "
  cold_pool -> eps_Pollock, 1, beta_cp_plk, 0,
  eps_Pollock -> eps_Pollock, 1, rho_plk, 0,
  cold_pool <-> cold_pool, 0, sigma_cp, 0.1,
  eps_Pollock <-> eps_Pollock, 0, sigma_plk, 0.5
"

## Initial run (to return simulate()) ----------------------

# Run with SEM
out0 <- ecostate(
  taxa = taxa, years = years, catch = catch, biomass = biomass, 
  PB = PB, QB = QB, B = B, DC = DC, EE = EE, X = X, 
  type = type, U = U, fit_Q = fit_Q, fit_B0 = fit_B0, fit_B = fit_B, fit_EE = fit_EE,
  sem = sem, covariates = covariates,
  settings = settings, control = control
)

## Loop over parameter grid, fit models -------------------

n_reps <- 3

beta_grid <- as.matrix(expand.grid(
  beta_cp_plk = seq(-0.5, 0.5, 0.25), 
  rho_plk = seq(-0.2, 0.8, 0.2)
))

par <- out0$internal$parhat
parhat <- array(NA, dim = c(nrow(beta_grid), ncol = length(out0$opt$par), n_reps))

set.seed(100)

for (i in 1:nrow(beta_grid)) {
  
  cat(paste0("beta ", i, "/", nrow(beta_grid), " ..........\n"))
  
  for (j in 1:n_reps) {
    
    cat(paste0("  rep ", j, "/", n_reps, "\n"))
    
    par$beta[1:2] <- beta_grid[i,] 
    sim_ij <- out0$simulator(par)
    
    catch <- na.omit(cbind(
      expand.grid( "Year" = rownames(sim_ij$Cobs_ti),
                   "Taxon" = colnames(sim_ij$Cobs_ti) ),
      "Mass" = as.vector(sim_ij$Cobs_ti)
    ))
    
    biomass <- na.omit(cbind(
      expand.grid( "Year" = rownames(sim_ij$Bobs_ti),
                   "Taxon" = colnames(sim_ij$Bobs_ti) ),
      "Mass" = as.vector(sim_ij$Bobs_ti)
    ))
    
    covariates <- sim_ij$covariates
    
    fit_ij <- ecostate(
      taxa = taxa, years = years, catch = catch, biomass = biomass, 
      PB = PB, QB = QB, B = B, DC = DC, EE = EE, X = X, 
      type = type, U = U, fit_Q = fit_Q, fit_B0 = fit_B0, fit_B = fit_B, fit_EE = fit_EE,
      sem = sem, covariates = covariates,
      settings = settings, control = control
    )
    
    parhat[i,,j] <- fit_ij$opt$par
    
  }
  
}

parhat_long <- reshape2::melt(parhat)
colnames(parhat_long) <- c("index", "parameter", "run", "estimate")
parhat_long$true <- ifelse(parhat_long$parameter == 5, beta_grid[parhat_long$index, 1], ifelse(parhat_long$parameter == 6, beta_grid[parhat_long$index,2], out0$opt$par[parhat_long$parameter]))

## Comparison plots ---------------------------------------

parhat_long |> 
  filter(parameter == 5) |> 
  ggplot(aes(true, estimate, group = true)) + 
  geom_boxplot(fill = "black", alpha = 0.25, color = "grey40") + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  ggtitle(expression(paste("cold pool " %->% " ", epsilon["Pollock"], ", lag 1")))

ggsave("cold_pool_sim.png", height = 4, width = 6, units = "in", dpi = 300)

parhat_long |> 
  filter(parameter == 6) |> 
  ggplot(aes(true, estimate, group = true)) + 
  geom_boxplot(fill = "black", alpha = 0.25, color = "grey40") + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  ggtitle(expression(paste(epsilon["Pollock"] %->% " ", epsilon["Pollock"], ", lag 1")))

ggsave("epsilon_pollock_AR1.png", height = 4, width = 6, units = "in", dpi = 300)

parhat_long |> 
  filter(parameter != 5 & parameter != 6) |> 
  ggplot(aes(factor(parameter), estimate)) + 
  geom_point(aes(y = true), color = "red", size = 4) + 
  geom_point(alpha = 0.25) + 
  ggtitle("other parameters (constant)")

ggsave("other_parameters.png", height = 4, width = 6, units = "in", dpi = 300)


# Simulation testing, epsilon + nu ------------------------

# DSEM paths
sem <- "
  eps_Pollock <-> eps_Pollock, 0, sigma_plk, 0.5,
  cold_pool -> nu_Arrowtooth, 0, beta_cp_atf, 0,
  eps_Pollock -> eps_Pollock, 1, rho_plk, 0,
  cold_pool <-> cold_pool, 0, sigma_cp, 0.1
"

## Initial run (to return simulate()) ---------------------

# Run with SEM
out0 <- ecostate(
  taxa = taxa, years = years, catch = catch, biomass = biomass, 
  PB = PB, QB = QB, B = B, DC = DC, EE = EE, X = X, 
  type = type, U = U, fit_Q = fit_Q, fit_B0 = fit_B0, fit_B = fit_B, fit_EE = fit_EE,
  sem = sem, covariates = covariates,
  settings = settings, control = control
)

## Loop over parameter grid, fit models -------------------

n_reps <- 3

beta_grid <- as.matrix(expand.grid(
  beta_cp_atf = seq(-0.5, 0.5, 0.25), 
  rho_plk = seq(-0.2, 0.8, 0.2)
))

par <- out0$internal$parhat
parhat <- array(NA, dim = c(nrow(beta_grid), ncol = length(out0$opt$par), n_reps))

set.seed(100)

for (i in 1:nrow(beta_grid)) {
  
  cat(paste0("beta ", i, "/", nrow(beta_grid), " ..........\n"))
  
  for (j in 1:n_reps) {
    
    cat(paste0("  rep ", j, "/", n_reps, "\n"))
    
    par$beta[2:3] <- beta_grid[i,] 
    sim_ij <- out0$simulator(par)
    
    catch <- na.omit(cbind(
      expand.grid( "Year" = rownames(sim_ij$Cobs_ti),
                   "Taxon" = colnames(sim_ij$Cobs_ti) ),
      "Mass" = as.vector(sim_ij$Cobs_ti)
    ))
    
    biomass <- na.omit(cbind(
      expand.grid( "Year" = rownames(sim_ij$Bobs_ti),
                   "Taxon" = colnames(sim_ij$Bobs_ti) ),
      "Mass" = as.vector(sim_ij$Bobs_ti)
    ))
    
    covariates <- sim_ij$covariates
    
    fit_ij <- ecostate(
      taxa = taxa, years = years, catch = catch, biomass = biomass, 
      PB = PB, QB = QB, B = B, DC = DC, EE = EE, X = X, 
      type = type, U = U, fit_Q = fit_Q, fit_B0 = fit_B0, fit_B = fit_B, fit_EE = fit_EE,
      sem = sem, covariates = covariates,
      settings = settings, control = control
    )
    
    parhat[i,,j] <- fit_ij$opt$par
    
  }
  
}

parhat_long <- reshape2::melt(parhat)
colnames(parhat_long) <- c("index", "parameter", "run", "estimate")
parhat_long$true <- ifelse(parhat_long$parameter == 6, beta_grid[parhat_long$index, 1], ifelse(parhat_long$parameter == 7, beta_grid[parhat_long$index,2], out0$opt$par[parhat_long$parameter]))


parhat_long |> 
  filter(parameter == 6) |> 
  ggplot(aes(true, estimate, group = true)) + 
  geom_boxplot(fill = "black", alpha = 0.25, color = "grey40") + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  ggtitle(expression(paste("cold pool " %->% " ", nu["Arrowtooth"], ", lag 0")))

ggsave("atf_nu_sim.png", height = 4, width = 6, units = "in", dpi = 300)

parhat_long |> 
  filter(parameter == 7) |> 
  ggplot(aes(true, estimate, group = true)) + 
  geom_boxplot(fill = "black", alpha = 0.25, color = "grey40") + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  ggtitle(expression(paste(epsilon["Pollock"] %->% " ", epsilon["Pollock"], ", lag 1")))

parhat_long |> 
  filter(parameter != 6 & parameter != 7) |> 
  ggplot(aes(factor(parameter), estimate)) + 
  geom_point(aes(y = true), color = "red", size = 4) + 
  geom_point(alpha = 0.25) + 
  ggtitle("other parameters (constant)")

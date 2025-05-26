
library("dsem")
library("RTMB")
library("ecostate")
library("tidyverse")

# Load data -----------------------------------------------

data(gulf_of_alaska)
goa <- gulf_of_alaska

goa$taxa <- gsub(" ", ".", goa$taxa)
names(goa$type) <- names(goa$B) <- names(goa$P_over_B) <- names(goa$Q_over_B) <- names(goa$EE) <- goa$taxa
colnames(goa$Diet_proportions) <- rownames(goa$Diet_proportions) <- goa$taxa
goa$catch_data$Taxon <- gsub(" ", ".", goa$catch_data$Taxon)
goa$biomass_data$Taxon <- gsub(" ", ".", goa$biomass_data$Taxon)

names(goa$stanza_groups) <- names(goa$Amax) <- gsub(" ", ".", names(goa$stanza_groups))
goa$stanza_groups <- gsub(" ", ".", goa$stanza_groups)
names(goa$K) <- names(goa$d) <- names(goa$Amat) <- names(goa$Wmat) <- names(goa$Wmatslope) <- names(goa$agecomp_data) <- gsub(" ", ".", names(goa$K))
names(goa$Leading) <- gsub(" ", ".", names(goa$Leading))

attach(goa)

# Inputs, Parameters to estimate --------------------------

# Biomass-dynamics steps
n_step = 50

# Age-structured dynamics steps
STEPS_PER_YEAR = 24

# Constant expected recruitment (matching assessment models)
# So h = 0.999 and back-calculate SpawnX
SpawnX = c( "Walleye.pollock" = 4 / (5-1/0.999), "Sablefish" = 4 / (5-1/0.999) )

# Default vulnerability as fixed or starting values
X = array( 2, dim = rep(length(taxa),2), dimnames = list(Prey=taxa,Predator=taxa) )

# Unassimilated food
U = array( 0.2, dim=length(taxa), dimnames=list(taxa) )

# Estimate catchability coefficient for all surveys
fit_Q = c( "Walleye.pollock.adult", "Sablefish.adult", "Euphausiids", "Large.copepods" )

# Fit equilibrium biomass for age-structured populations
# using fishery as depletion experiment to identify scale
fit_B = c("Walleye.pollock.juv", "Sablefish.adult")

# Fit PB with a prior
fit_PB = c( "Walleye.pollock.adult", "Sablefish.adult" )

# Fit vulnerability for adult age-structured populations with a prior
fit_X = c("Walleye pollock adult", "Sablefish adult")

# SEM - biomass and recruitment deviations
sem = "
  eps_Euphausiids <-> eps_Euphausiids, 0, tau_Eup, 1,
  eps_Large.copepods <-> eps_Large.copepods, 0, tau_Cop, 1, 
  phi_Walleye.pollock <-> phi_Walleye.pollock, 0, psi_Plk, 1,
  phi_Sablefish <-> phi_Sablefish, 0, psi_Sab, 1
"

# Priors
log_prior = list(
  
  # Normal(0,0.5) prior on diag(Xprime_ij), where X = exp(Xprime) + 1
  #diag(Xprime_ij) ~ dnorm(mean = 0, sd = 0.5),
  
  # Normal prior on log(q) for adult pollock, matching stock assessment
  logq_i['Walleye.pollock.adult'] ~ dnorm(mean = log(0.85), sd = 0.1),
  
  # Tight normal prior on log(PB) for sablefish, matching assessment M value
  logPB_i["Sablefish.adult"] ~ dnorm(mean = log(0.1), sd = 0.1),
  
  # Tight normal prior on log(PB) for pollock, matching assessment M value
  logPB_i['Walleye.pollock.adult'] ~ dnorm(mean = log(0.3), sd = 0.1)#, 
  
  #abs(c(tau_Eup, tau_Cop, psi_Plk, psi_Sab)) ~ dexp(rate = 2)
  
)

# SEM ML epsilon + phi implementation ---------------------

control <- ecostate_control( 
  n_steps = n_step,   
  profile = NULL,
  random = NULL,
  derived_quantities = c(), 
  getsd = FALSE,
  silent = FALSE,
  verbose = TRUE, 
  trace = TRUE
)

settings_dgmrf <- stanza_settings(
  taxa = taxa, stanza_groups = stanza_groups,
  K = K, Wmat = Wmat, Amat = Amat, d = d, Amax = Amax, SpawnX = SpawnX,
  STEPS_PER_YEAR = STEPS_PER_YEAR, comp_weight = "multinom", Leading = Leading,
  Wmatslope = Wmatslope
)

dgmrf_fit <- ecostate(
  taxa = taxa, years = years, type = type,
  catch = catch_data, biomass = biomass_data, agecomp = agecomp_data,
  PB = P_over_B, QB = Q_over_B, DC = Diet_proportions, B = B, EE = EE, X = X, U = U,
  fit_B = fit_B, fit_Q = fit_Q, fit_PB = fit_PB, sem = sem,
  log_prior = log_prior, settings = settings_dgmrf, control = control
)

# Non-SEM ML implementation -------------------------------

# Fit biomass-dynamics process errors for zooplankton
fit_eps = c( "Euphausiids", "Large.copepods" )

# Fit recruitment deviations for age-structured populations
fit_phi = c("Sablefish", "Walleye.pollock")

settings_dnorm <- stanza_settings(
  taxa = taxa, stanza_groups = stanza_groups,
  K = K, Wmat = Wmat, Amat = Amat, d = d, Amax = Amax, SpawnX = SpawnX,
  STEPS_PER_YEAR = STEPS_PER_YEAR, comp_weight = "multinom", Leading = Leading,
  Wmatslope = Wmatslope, fit_phi = fit_phi
)

dnorm_fit <- ecostate(
  taxa = taxa, years = years, type = type,
  catch = catch_data, biomass = biomass_data, agecomp = agecomp_data,
  PB = P_over_B, QB = Q_over_B, DC = Diet_proportions, B = B, EE = EE, X = X, U = U,
  fit_B = fit_B, fit_Q = fit_Q, fit_PB = fit_PB, fit_eps = fit_eps,
  log_prior = log_prior, settings = settings_dnorm, control = control
)

# Fixed run -----------------------------------------------

# Fix SD of biomass and recruitment errors (penalized likelihood)
sem = "
  eps_Euphausiids <-> eps_Euphausiids, 0, NA, 1,
  eps_Large.copepods <-> eps_Large.copepods, 0, NA, 1, 
  phi_Walleye.pollock <-> phi_Walleye.pollock, 0, NA, 1,
  phi_Sablefish <-> phi_Sablefish, 0, NA, 1
"

dgmrf_penalized <- ecostate(
  taxa = taxa, years = years, type = type,
  catch = catch_data, biomass = biomass_data, agecomp = agecomp_data,
  PB = P_over_B, QB = Q_over_B, DC = Diet_proportions, B = B, EE = EE, X = X, U = U,
  fit_B = fit_B, fit_Q = fit_Q, fit_PB = fit_PB, sem = sem,
  log_prior = log_prior, settings = settings_dgmrf, control = control
)

# Fixed run, non-SEM implementation -----------------------

control$nlminb_loops <- 0

dnorm_penalized0 <- ecostate(
  taxa = taxa, years = years, type = type,
  catch = catch_data, biomass = biomass_data, agecomp = agecomp_data,
  PB = P_over_B, QB = Q_over_B, DC = Diet_proportions, B = B, EE = EE, X = X, U = U,
  fit_B = fit_B, fit_Q = fit_Q, fit_PB = fit_PB, fit_eps = fit_eps,
  log_prior = log_prior, settings = settings_dnorm, control = control
)

map <- dnorm_penalized0$tmb_inputs$map
tmb_par <- dnorm_penalized0$obj$env$parList()

# Penalized likelihood:  Fixed SD for process errors = 1
map$logtau_i = factor(rep(NA,length(map$logtau_i)))
tmb_par$logtau_i = ifelse( is.na(tmb_par$logtau_i), NA, log(1) ) 

# Penalized likelihood:  Fixed SD for recruitment deviations = 1
map$logpsi_g2 = factor(rep(NA,length(map$logpsi_g2)))
tmb_par$logpsi_g2 = ifelse( is.na(tmb_par$logpsi_g2), NA, log(1) ) 

pcontrol <- control
pcontrol$nlminb_loops <- 1
pcontrol$map <- map
pcontrol$tmb_par <- tmb_par

dnorm_penalized <- ecostate(
  taxa = taxa, years = years, type = type,
  catch = catch_data, biomass = biomass_data, agecomp = agecomp_data,
  PB = P_over_B, QB = Q_over_B, DC = Diet_proportions, B = B, EE = EE, X = X, U = U,
  fit_B = fit_B, fit_Q = fit_Q, fit_PB = fit_PB, fit_eps = fit_eps,
  log_prior = log_prior, settings = settings_dnorm, control = pcontrol
)

# Comparison ----------------------------------------------

## Likelihood, same parameters ----------------------------

dgmrf_par <- dgmrf_fit$opt$par
dgmrf_par[which(names(dgmrf_par) == "beta")] <- 1

dnorm_par <- dnorm_fit$opt$par
dnorm_par[!(names(dnorm_par) %in% c("logtau_i", "logpsi_g2"))] <- dgmrf_par[names(dgmrf_par) != "beta"]
dnorm_par[names(dnorm_par) %in% c("logtau_i", "logpsi_g2")] <- log(1)

cor(dgmrf_par[!(names(dgmrf_par) == "beta")], dnorm_par[!(names(dnorm_par) %in% c("logtau_i", "logpsi_g2"))])

dgmrf_fit$obj$fn(dgmrf_par) == dnorm_fit$obj$fn(dnorm_par)

## Penalized models ---------------------------------------

## Penalized parameters match 1-to-1
cor(dnorm_penalized$opt$par, dgmrf_penalized$opt$par)

## Accurate to 2 decimal places
all(round(dnorm_penalized$opt$par, 2) == round(dgmrf_penalized$opt$par, 2))

## Plot estimated process errors
plot(dnorm_penalized$internal$parhat$epsilon_ti[,7], dgmrf_penalized$internal$parhat$epsilon_ti[,7]); abline(0, 1)
plot(dnorm_penalized$internal$parhat$epsilon_ti[,11], dgmrf_penalized$internal$parhat$epsilon_ti[,11]); abline(0, 1)
plot(dnorm_penalized$internal$parhat$phi_tg2[,1], dgmrf_penalized$internal$parhat$phi_tg2[,1]); abline(0, 1)
plot(dnorm_penalized$internal$parhat$phi_tg2[,2], dgmrf_penalized$internal$parhat$phi_tg2[,2]); abline(0, 1)

## Non-penalized models -----------------------------------

## Non-penalized - standard deviations do not match
abs(dgmrf_fit$internal$parhat$beta)
exp(c(dnorm_fit$internal$parhat$logtau_i[c(7, 11)], dnorm_fit$internal$parhat$logpsi_g2))

## Process errors - no agreement
cor(dnorm_fit$internal$parhat$epsilon_ti[,7], dgmrf_fit$internal$parhat$epsilon_ti[,7]) # 0.148
cor(dnorm_fit$internal$parhat$epsilon_ti[,11], dgmrf_fit$internal$parhat$epsilon_ti[,11]) # -0.099
cor(dnorm_fit$internal$parhat$phi_tg2[,1], dgmrf_fit$internal$parhat$phi_tg2[,1]) # 0.347
cor(dnorm_fit$internal$parhat$phi_tg2[,2], dgmrf_fit$internal$parhat$phi_tg2[,2]) # 0.012

## Other parameters - pretty good agreement (0.994)
cor (
  dgmrf_fit$opt$par[!(names(dgmrf_fit$opt$par) %in% c("epsilon_ti", "phi_tg2", "beta", "logtau_i", "logpsi_g2"))], 
  dnorm_fit$opt$par[!(names(dgmrf_fit$opt$par) %in% c("epsilon_ti", "phi_tg2", "beta", "logtau_i", "logpsi_g2"))]
)

savedir <- "C:/Users/goodm/Google Drive/UW Res Sci/EcoState/scratch/"
saveRDS(dgmrf_fit, paste0(savedir, "dgmrf_fit.rds"))
saveRDS(dnorm_fit, paste0(savedir, "dnorm_fit.rds"))
saveRDS(dgmrf_penalized, paste0(savedir, "dgmrf_penalized.rds"))
saveRDS(dnorm_penalized, paste0(savedir, "dnorm_penalized.rds"))

# Fit penalized model with temperature effect ---------------

## SST summarized from monthly time series https://shinyfin.psmfc.org/ak-sst-mhw/
goa_sst <- structure(list(
  year = c(
    1985, 1986, 1987, 1988, 1989, 1990, 1991, 
    1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 
    2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 
    2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024
  ), 
  sst = c(
    7.232, 7.699, 7.997, 7.661, 7.571, 7.742, 7.356, 7.625, 
    8.056, 7.685, 7.595, 7.665, 8.284, 8.174, 7.323, 7.733, 7.985, 
    7.862, 8.605, 8.369, 8.577, 7.738, 7.437, 7.105, 7.416, 7.819, 
    7.534, 7.05, 7.821, 8.529, 8.793, 9.291, 8.328, 8.411, 9.033, 
    8.299, 7.901, 8.074, 7.886, 8.043
  )), 
  row.names = c(NA, -40L), class = "data.frame"
)

sem = "
  eps_Euphausiids <-> eps_Euphausiids, 0, NA, 1,
  eps_Large.copepods <-> eps_Large.copepods, 0, NA, 1, 
  phi_Walleye.pollock <-> phi_Walleye.pollock, 0, NA, 1,
  phi_Sablefish <-> phi_Sablefish, 0, NA, 1, 
  sst -> phi_Walleye.pollock, 0, beta_sst_plk, 0
"

covariates <- matrix(NA, nrow = length(years), ncol = 1)
colnames(covariates) <- "sst"
covariates[match(goa_sst$year, years)] <- goa_sst$sst

sst_ctrl <- ecostate_control(
  n_steps = n_step, profile = NULL, random = c("covariates"), getsd = TRUE, 
  silent = FALSE, verbose = TRUE, trace = TRUE
)

sst_fit <- ecostate(
  taxa = taxa, years = years, type = type,
  catch = catch_data, biomass = biomass_data, agecomp = agecomp_data,
  PB = P_over_B, QB = Q_over_B, DC = Diet_proportions, B = B, EE = EE, X = X, U = U,
  fit_B = fit_B, fit_Q = fit_Q, fit_PB = fit_PB, sem = sem, covariates = covariates,
  log_prior = log_prior, settings = settings_dgmrf, control = sst_ctrl
)

plot()
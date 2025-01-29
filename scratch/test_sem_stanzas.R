
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

# Fit biomass-dynamics process errors for zooplankton
fit_eps = c( "Euphausiids", "Large.copepods" )

# Fit recruitment deviations for age-structured populations
fit_phi = c("Walleye.pollock", "Sablefish")

# Fit equilibrium biomass for age-structured populations
# using fishery as depletion experiment to identify scale
fit_B = c("Walleye.pollock.juv", "Sablefish.adult")

# Fit PB with a prior
fit_PB = c( "Walleye.pollock.adult", "Sablefish.adult" )

# SEM - biomass and recruitment deviations
# Fix SD of biomass and recruitment errors (penalized likelihood)
sem = "
  eps_Euphausiids <-> eps_Euphausiids, 0, NA, 1,
  eps_Large.copepods <-> eps_Large.copepods, 0, NA, 1, 
  phi_Walleye.pollock <-> phi_Walleye.pollock, 0, NA, 1,
  phi_Sablefish <-> phi_Sablefish, 0, sigma, 1
"

# Priors
log_prior = list(
  
  # Normal prior on log(q) for adult pollock, matching stock assessment
  logq_i['Walleye.pollock.adult'] ~ dnorm(mean = log(0.85), sd = 0.1),
  
  # Tight normal prior on log(PB) for sablefish, matching assessment M value
  logPB_i["Sablefish.adult"] ~ dnorm(mean = log(0.1), sd = 0.1),
  
  # Tight normal prior on log(PB) for pollock, matching assessment M value
  logPB_i['Walleye.pollock.adult'] ~ dnorm(mean = log(0.3), sd = 0.1)
  
)

# Initial run ---------------------------------------------

settings <- stanza_settings(
  taxa = taxa, stanza_groups = stanza_groups,
  K = K, Wmat = Wmat, Amat = Amat, d = d, Amax = Amax, SpawnX = SpawnX,
  STEPS_PER_YEAR = STEPS_PER_YEAR, comp_weight = "multinom", Leading = Leading,
  Wmatslope = Wmatslope
)

control <- ecostate_control( 
  n_steps = n_step,   
  profile = NULL, # Penalized likelihood so use empty set
  random = NULL # Penalized likelihood so use empty set
)

out0 <- ecostate(
  taxa = taxa, years = years, type = type,
  catch = catch_data, biomass = biomass_data, agecomp = agecomp_data,
  PB = P_over_B, QB = Q_over_B, DC = Diet_proportions, B = B, EE = EE, X = X, U = U,
  fit_B = fit_B, fit_Q = fit_Q, fit_PB = fit_PB, sem = sem,
  log_prior = log_prior, settings = settings, # control = control,
  debug = TRUE
)

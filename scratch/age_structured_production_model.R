
library(ecostate)
library(JABBA)
library(spict)
set.seed(101)

# Time-interval
years = 1981:2020
n_years = length(years)

# Biology
r = 0.2
#MSY = 100
K = 1000
sigmaB = 0.5
B0 = K * exp(sigmaB*rnorm(1))
prod_func = c("Schaefer", "Fox")[2]

# Effort dynamics
Bequil = 0.4 * K
Brate = 0.2
sigmaE = 0.1
E0 = 0.01

# Survey
q = 1
sigmaQ = 0.1

#
Pbar_t = P_t = C_t = E_t = B_t = rep(NA, n_years)
B_t[1] = B0
E_t[1] = E0

#
for( t in 2:n_years ){
  if(prod_func=="Scaefer") Pbar_t[t] = B_t[t-1] + r * B_t[t-1] * (1 - B_t[t-1]/K)
  if(prod_func=="Fox") Pbar_t[t] = r * B_t[t-1] * log(K / B_t[t-1])
  P_t[t] = Pbar_t[t] * exp(sigmaB*rnorm(1))
  B_t[t] = B_t[t-1] + P_t[t]
  E_t[t] = E_t[t-1] * (B_t[t-1]/Bequil)^Brate * exp(sigmaE*rnorm(1))
  C_t[t] = B_t[t] * (1 - exp(-E_t[t]))
  B_t[t] = B_t[t] - C_t[t]
}
Bobs_t = q * B_t * exp(sigmaQ * rnorm(n_years))

# 
matplot( x=years, y=cbind(B_t,Bobs_t,C_t), type="l", log="y", lty="solid")

# Name taxa (optional, for illustration)
taxa = "target"
n_taxa = length(taxa)

# Ecopath-with-EcoSim parameters
# Diet matrix
DC_ij = array( 0, dim=c(1,1) )
PB_i = 0.1
QB_i = NA
EE_i = 1
B_i = 1
U_i = 0.2
type_i = "auto"
X_ij = array( 2, dim=c(1,1) )

# reformat to longform data-frame
Catch = na.omit(data.frame( "Mass" = C_t, "Year" = years, "Taxon" = taxa ))
Biomass = data.frame( "Mass" = Bobs_t, "Year" = years, "Taxon" = taxa )

# Inputs
taxa = c( "prey", "stage1", "stage2" )
n_taxa = length(taxa)
PB_i = c( 5, 1, 0.1 )
QB_i = c( NA, NA, 0.5 )
DC_ij = cbind( "prey"=c(0,0,0), "stage1"=c(1,0,0), "stage2"=c(1,0,0)  )
X_ij = matrix( 2, nrow=length(taxa), ncol=length(taxa) )
U_i = c( 0.2, 0.2, 0.2 )
type = c( "auto", "hetero", "hetero" )
EE_i = c( 1, NA, NA ) 
B_i = c( NA, NA, 1 )

# Settings: specify what parameters to estimate
fit_Q = "stage2"       # catchability coefficient
fit_B0 = c()      # non-equilibrium initial condition
fit_B = "stage2"       # equilibrium biomass

# Label EwE inputs for each taxon as expected (so users can easily change taxa)
names(PB_i) = names(QB_i) = names(B_i) = names(EE_i) = names(type) = names(U_i) = taxa
dimnames(DC_ij) = dimnames(X_ij) = list("Prey"=taxa, "Predator"=taxa)

# Fit using SEM -------------------------------------------

settings <- stanza_settings(
  taxa = taxa,
  stanza_groups = c("stage1"="age_structured", "stage2"="age_structured"),
  K = c("age_structured" = 0.1),
  Wmat = c("age_structured" = (2/3)^3),
  d = c("age_structured" = 2/3),
  Amax = c("stage1" = 2, "stage2" = 30),
  SpawnX = c("age_structured" = 2)
)

sem <- "
  eps_stage1 <-> eps_stage1, 0, tau1, 0.001, 
  eps_stage2 <-> eps_stage2, 0, tau2, 0.001
"

# Prior on process-error log-SD to stabilize model
log_prior <- list(
  log(c(tau1, tau2)) ~ dnorm(mean = log(0.2), sd = 1)
)

# Run model
sem_fit = ecostate( 
  taxa = taxa,
  years = years,
  catch = data.frame("Mass"=Catch[,1],"Year"=Catch[,2],"Taxon"="stage2"),
  biomass = data.frame("Mass"=Biomass[,1],"Year"=Biomass[,2],"Taxon"="stage2"),
  PB = PB_i,
  QB = QB_i,
  DC = DC_ij,
  B = B_i,
  EE = EE_i,
  X = X_ij,
  type = type,
  sem = sem,
  U = U_i,
  fit_B = fit_B,
  fit_Q = fit_Q,
  fit_B0 = fit_B0,
  log_prior = log_prior,
  control = ecostate_control(),
  settings = settings
)

# All parameters
sem_fit$opt$par

# SD of biomass deviations
abs(sem_fit$sem$Estimate)

# Fit using fit_eps ---------------------------------------

fit_eps = c("stage1", "stage2")

# Define priors
log_prior = list(
  logtau_i[c("stage1", "stage2")] ~ dnorm(mean = log(0.2), sd = 1)
)

# Run model
out_fiteps = ecostate( 
  taxa = taxa,
  years = years,
  catch = data.frame("Mass"=Catch[,1],"Year"=Catch[,2],"Taxon"="stage2"),
  biomass = data.frame("Mass"=Biomass[,1],"Year"=Biomass[,2],"Taxon"="stage2"),
  PB = PB_i,
  QB = QB_i,
  DC = DC_ij,
  B = B_i,
  EE = EE_i,
  X = X_ij,
  type = type,
  U = U_i,
  fit_B = fit_B,
  fit_Q = fit_Q,
  fit_eps = fit_eps,
  fit_B0 = fit_B0,
  log_prior = log_prior,
  control = ecostate_control(),
  settings = settings
)

out_fiteps$opt$par

exp(out_fiteps$opt$par[names(out_fiteps$opt$par) == "logtau_i"])


# devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\ecostate)', force=TRUE )
library(ecostate)

# load data
data(eastern_bering_sea)

# Reformat inputs
years = 1982:2021 # Catch only goes through 2021, and starting pre-data in 1982 doesn't play well with fit_B0
taxa = c( "Pollock", "Cod", "Arrowtooth", "Copepod", "Other_zoop", "Chloro", "NFS", "Krill", "Benthic_invert", "Benthos", "Detritus" )

# Define types
type_i = sapply( taxa, FUN=switch, "Detritus" = "detritus",
                                   "Chloro" = "auto",
                                   "hetero" )

# Starting values
U_i = EE_i = B_i = array( NA, dim=length(taxa),
                    dimnames=list(names(eastern_bering_sea$P_over_B)))
B_i[c("Pollock","Cod", "Arrowtooth", "NFS")] = c(1, 1, 0.5, 0.02)
EE_i[] = 1
U_i[] = 0.2

# Define default vulnerability, except for primary producers
X_ij = array( 2, dim=c(length(taxa),length(taxa)) )
dimnames(X_ij) = list(names(B_i),names(B_i))
X_ij[,'Chloro'] = 91

# Define parameters to estimate
fit_Q = c("Pollock", "Copepod", "Chloro", "Other_zoop", "Krill")
fit_B0 = c("Pollock", "Cod", "Arrowtooth", "NFS")
fit_B = c("Cod", "Arrowtooth", "NFS", "Pollock")

# Define process errors to estimate
# Only estimating Pollock to speed up demonstration
fit_eps = "Pollock"
#fit_eps = c()

# Which taxa to include
taxa_to_include = c( "NFS", "Pollock", "Copepod", "Chloro", "Krill" )
# To run full model use:
# taxa_to_include = taxa

# Define priors
log_prior = function(p){
  # Prior on process-error log-SD to stabilize model
  logp = sum(dnorm( p$logtau_i, mean=log(0.2), sd=1, log=TRUE ), na.rm=TRUE)
}

# Run model
out0 = ecostate( taxa = taxa_to_include,
                years = years,
                catch = eastern_bering_sea$Catch,
                biomass = eastern_bering_sea$Survey,
                PB = eastern_bering_sea$P_over_B,
                QB = eastern_bering_sea$Q_over_B,
                DC = eastern_bering_sea$Diet_proportions,
                B = B_i,
                EE = EE_i,
                U = U_i,
                type = type_i,
                X = X_ij,
                fit_B = fit_B,
                fit_Q = fit_Q,
                fit_eps = fit_eps,
                fit_B0 = fit_B0,
                log_prior = log_prior,
                control = ecostate_control( n_steps = 20, # using 15 by default
                                            start_tau = 0.2,
                                            tmbad.sparse_hessian_compress = 0,
                                            getsd = TRUE,
                                            #eval.max = 10,
                                            #iter.max = 10,
                                            trace = 1 ))

# print output
out0

#
set.seed(101)
sim = out0$simulator()
matplot( x = out0$internal$years,
         y = sim$B_ti,
         type = "l",
         lwd = 2,
         log = "y",
         col = rainbow(length(out$internal$taxa)) )
legend( "left",
        fill = rainbow(length(out$internal$taxa)),
        legend=out$internal$taxa )

# Run model
catch = na.omit(cbind(
  expand.grid( "Year" = rownames(sim$Cobs_ti),
               "Taxon" = colnames(sim$Cobs_ti) ),
  "Mass" = as.vector(sim$Cobs_ti)
))
biomass = na.omit(cbind(
  expand.grid( "Year" = rownames(sim$Bobs_ti),
               "Taxon" = colnames(sim$Bobs_ti) ),
  "Mass" = as.vector(sim$Bobs_ti)
))
out = ecostate( taxa = taxa_to_include,
                years = years,
                catch = catch,
                biomass = biomass,
                PB = eastern_bering_sea$P_over_B,
                QB = eastern_bering_sea$Q_over_B,
                DC = eastern_bering_sea$Diet_proportions,
                B = B_i,
                EE = EE_i,
                U = U_i,
                type = type_i,
                X = X_ij,
                fit_B = fit_B,
                fit_Q = fit_Q,
                fit_eps = fit_eps,
                fit_B0 = fit_B0,
                log_prior = log_prior,
                control = ecostate_control( n_steps = 20, # using 15 by default
                                            start_tau = 0.2,
                                            tmbad.sparse_hessian_compress = 0,
                                            getsd = TRUE,
                                            #eval.max = 10,
                                            #iter.max = 10,
                                            trace = 1 ))


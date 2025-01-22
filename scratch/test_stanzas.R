
library("checkmate")
library("Matrix")
library("dsem")
library("RTMB")
library("ecostate")
library("tidyverse")

# Load data -----------------------------------------------

data( whitehouse_2021 )
Ecopath = whitehouse_2021$Ecopath
Diet = whitehouse_2021$Diet
stanzas = whitehouse_2021$stanzas
stanza_groups = whitehouse_2021$stanza_groups

data( eastern_bering_sea )
Catch = eastern_bering_sea$Catch
catch = data.frame( "Year"=Catch$Year, "Mass"=Catch$Mass, "Taxon"=Catch$Taxon )
catch$Taxon = sapply( catch$Taxon, FUN=switch, "Arrowtooth" = "Arrowtooth flounder adult",
                      "Cod" = "Pacific cod adult",
                      "Pollock" = "Walleye pollock adult" )

#
Fish = eastern_bering_sea$Survey
biomass = data.frame(
  "Mass" = Fish$Mass,
  "Year" = Fish$Year,
  "Taxon" = sapply( Fish$Taxon, FUN=switch, "Arrowtooth" = "Arrowtooth flounder adult",
                    "Cod" = "Pacific cod adult",
                    "Pollock" = "Walleye pollock adult",
                    NA)
)
biomass = na.omit(biomass)

Ecopath$Functional.group <- gsub(" ", ".", Ecopath$Functional.group)
which_ecopath = setdiff(1:nrow(Ecopath), 66)
taxa = c( Ecopath[which_ecopath,'Functional.group'], "benthic_detritus" )
type = sapply( taxa, FUN = switch, 
               "benthic_detritus" = "detritus", 
               "Large.phytoplankton" = "auto", 
               "Small.phytoplankton" = "auto", 
               "hetero" )
PB = c( Ecopath[which_ecopath,'P.B'], 0.5 )
QB = c( Ecopath[which_ecopath,'Q.B'], NA )
B = c( Ecopath[which_ecopath,'Biomass'], NA )
EE = c( Ecopath[which_ecopath,'EE'], 0.5 )
U = c( rep(0.2,length(which_ecopath)), 0 )

# Modify EE
EE = ifelse( taxa %in% c("benthic_detritus","Walleye.pollock.adult"), EE, NA )

# Diet compositions - subset to included species
which_cols = which_ecopath + 1
which_rows = c( which_ecopath, 72 )
DC = cbind( Diet[which_rows, which_cols], "benthic_detritus"=0 )
DC = ifelse( is.na(DC), 0, as.matrix(DC) )
X = array(2, dim=rep(length(taxa),2) )

# Stanza groups
stgroups = stanza_groups[,'StanzaGroup'][stanzas[,'StGroupNum']]
Amax = ceiling((stanzas[,'Last']+1) / 12)

# Life history parameters
K = stanza_groups[,'VBGF_Ksp']
d = stanza_groups[,'VBGF_d']
Wmat = stanza_groups[,'Wmat']
SpawnX = rep(2, length(K))

# Add names
names(PB) = names(QB) = names(B) = names(EE) = names(U) = names(type) = taxa
dimnames(DC) = dimnames(X) = list( taxa, taxa )
names(Amax) = names(stgroups) = taxa[stanzas[,'GroupNum']]
names(K) = names(d) = names(Wmat) = names(SpawnX) = unique(stgroups)

# Eliminate juvenile values
which_juv = stanzas[which(!stanzas[,'Leading']),'GroupNum']
QB[which_juv] = B[which_juv] = NA

# Remove spaces for dsem
catch$Taxon <- gsub(" ", ".", catch$Taxon)
biomass$Taxon <- gsub(" ", ".", biomass$Taxon)

# Age comps from survey data
load("C:/Users/goodm/Documents/R Projects/EBS_range_projections/data/trawl_surveys_size_binned/ebs.srvy98.plk.cpue_data.Rdata")
pollock <- cpue_data
load("C:/Users/goodm/Documents/R Projects/EBS_range_projections/data/trawl_surveys_size_binned/ebs.srvy98.atf.cpue_data.Rdata")
flounder <- cpue_data

vbgf_age <- function(Lt, Linf = 91, K = 0.192) -log((Linf - Lt)/Linf)/K

pollock <- pollock$CPUE_station_bin_yr |> 
  filter(BIN < 900) |> 
  mutate(age = floor(vbgf_age(BIN/10))) |> 
  group_by(YEAR, age) |> 
  summarize(cpue = sum(haul_CPUE_NUMKM2), .groups = "drop") |> 
  group_by(YEAR) |> 
  mutate(prop = round(cpue/sum(cpue) * 100)) |> 
  mutate(age = factor(age, levels = 0:34)) |> 
  ungroup() |> 
  complete(YEAR, age, fill = list(prop = 0)) |> 
  select(-cpue) |>
  pivot_wider(names_from = "age", values_from = "prop") |> 
  column_to_rownames("YEAR") |> 
  as.matrix()

flounder <- flounder$CPUE_station_bin_yr |> 
  filter(BIN < 840) |> 
  mutate(age = floor(vbgf_age(BIN/10, Linf = 84, K = 0.08775))) |> 
  group_by(YEAR, age) |> 
  summarize(cpue = sum(haul_CPUE_NUMKM2), .groups = "drop") |> 
  group_by(YEAR) |> 
  mutate(prop = round(cpue/sum(cpue) * 100)) |> 
  filter(age <= 34) |> 
  mutate(age = factor(age, levels = 0:34)) |> 
  ungroup() |> 
  complete(YEAR, age, fill = list(prop = 0)) |> 
  select(-cpue) |>
  pivot_wider(names_from = "age", values_from = "prop") |> 
  column_to_rownames("YEAR") |> 
  as.matrix()
  
agecomp <- list("W.Pollock" = pollock, "Arrowtooth" = flounder)

# Parameters to estimate --------------------------------------------

# Fit catchability coefficient for adult pollock
fit_Q = c( "Walleye.pollock.adult" )

# fit equilibrium biomass for adult arrowtooth and cod
fit_B = c( "Arrowtooth.flounder.adult", "Pacific.cod.adult" )

sem = "
  cold_pool -> phi_W.Pollock, 1, beta_cp_plk, 0,
  phi_W.Pollock -> phi_W.Pollock, 1, rho_plk, 0,
  cold_pool <-> cold_pool, 0, sigma_cp, 0.1
  phi_W.Pollock <-> phi_W.Pollock, 0, sigma_phi_plk, 0.1
"

# DSEM additional covariates
covariates <- structure(
  c(0.349, 0.278, 0.481, 0.491, 0.462, 0.3, 0.381, 0.181, 
    0.298, 0.379, 0.57, 0.223, 0.379, 0.59, 0.373, 0.56, 0.27, 0.708, 
    0.515, 0.424, 0.379, 0.292, 0.379, 0.351, 0.602, 0.682, 0.69, 
    0.7, 0.702, 0.544, 0.722, 0.615, 0.335, 0.355, 0.27, 0.475, 0.055, 
    0.187, 0.535, 0.295), 
  dim = c(40L, 1L), 
  dimnames = list(NULL, "cold_pool")
)

# Run model -----------------------------------------------

settings <- stanza_settings(
  taxa = taxa, stanza_groups = stgroups,
  K = K, Wmat = Wmat, d = d, Amax = Amax, SpawnX = SpawnX,
  STEPS_PER_YEAR = 12
)

out = ecostate(
  taxa = taxa, years = 1981:2020,
  catch = catch, biomass = biomass,
  PB = PB, QB = QB, DC = DC, B = B, EE = EE, X = X,
  type = type, U = U, fit_B = fit_B, fit_Q = fit_Q,
  sem = sem, covariates = covariates,
  settings = settings, 
  debug = FALSE
)

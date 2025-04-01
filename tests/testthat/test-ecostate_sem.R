
# EBS test inputs -----------------------------------------

data("eastern_bering_sea", package = "ecostate")
ebs <- eastern_bering_sea

# Data
taxa <- names(ebs$P_over_B)
years <- sort(unique(ebs$Survey$Year))
type <- sapply(taxa, switch, "Detritus" = "detritus", "Chloro" = "auto", "hetero")

# Initial equilibrium biomass and ecotrophic efficiencies
EE <- B <- U <- setNames(rep(NA, length(taxa)), taxa)
B[] <- aggregate(ebs$Survey$Mass, list(Taxon = factor(ebs$Survey$Taxon, levels = taxa)), mean, drop = FALSE)[,2]
EE[] <- 1
U[] <- 0.2

# Default vulnerability, except producers
X <- array(2, dim = c(length(taxa), length(taxa)))
dimnames(X) <- list(taxa, taxa)
X[,'Chloro'] <- 91

# Taxa to include
taxa <- c("Pollock", "Krill", "Arrowtooth", "Cod")

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

# Tests ---------------------------------------------------

test_that(
  "i.i.d. dgmrf eps/nu process errors match dnorm", {
    
    sem <- "
      eps_Pollock <-> eps_Pollock, 0, tau_plk, 0.5, 
      nu_Arrowtooth <-> nu_Arrowtooth, 0, sigma_atf, 0.5
    "
    
    # Base implementation
    dnorm_mod <- ecostate(
      taxa = taxa, years = years, catch = ebs$Catch, biomass = ebs$Survey, 
      PB = ebs$P_over_B, QB = ebs$Q_over_B, B = B, DC = ebs$Diet_proportions, 
      EE = EE, X = X, type = type, U = U, fit_B0 = "Pollock",
      control = ecostate_control(nlminb_loops = 0, getsd = FALSE, start_tau = 0.5), 
      fit_eps = "Pollock", fit_nu = "Arrowtooth"
    )
    
    # SEM implementation
    dgmrf_mod <- ecostate(
      taxa = taxa, years = years, catch = ebs$Catch, biomass = ebs$Survey, 
      PB = ebs$P_over_B, QB = ebs$Q_over_B, B = B, DC = ebs$Diet_proportions, 
      EE = EE, X = X, type = type, U = U, fit_B0 = "Pollock",
      control = ecostate_control(nlminb_loops = 0, getsd = FALSE), 
      sem = sem
    )
    
    # Store last best pars back to obj$env to speed up likelihood evaluation
    
    epsilon_ti <- c(
      0, -4.36393, 9.3658, 6.26889, 6.85014, 6.71161, 6.66377, 6.57245, 
      6.72672, 6.5751, 6.78338, 6.70492, 6.72848, 6.66424, 6.53801, 
      6.62825, 6.52334, 6.49804, 6.59593, 6.82203, 6.9659, 6.96004, 
      6.8289, 7.00278, 6.83973, 6.83431, 6.1379, 5.84587, 5.76003, 
      5.93233, 5.96902, 5.85114, 5.74212, 5.5586, 3.14531, 5.86583, 
      4.8834, 5.12842, 1.86746, 9.06086, 0.61628, 9.55358
    )
    
    nu_ti <- c(
      0, -2.70888, -3.3762, -3.30752, -3.25986, -3.20315, -3.1454, 
      -3.08645, -3.02677, -2.96364, -2.89834, -2.82886, -2.75755, -2.68375, 
      -2.61152, -2.53677, -2.45821, -2.37788, -2.28919, -2.20271, -2.11758, 
      -2.02197, -1.93165, -1.78828, -1.69276, -1.53434, -1.35108, -1.17176, 
      -0.93038, -0.67238, -0.59391, -0.38448, -0.35125, 5.47219, -1.03528, 
      -0.91892, -0.79909, 4.46919, -1.0059, -1.04759, -0.79686, -1.4701
    )
    
    logF_ti <- c(
      1.72168, 1.94369, 1.97134, 1.97808, 1.9822, 1.97224, 1.97503, 
      1.98864, 1.98469, 2.00688, 1.99043, 2.00394, 1.98451, 1.99348, 
      1.97522, 1.98626, 1.95806, 1.97413, 2.01757, 2.03106, 2.01633, 
      2.03258, 2.01741, 2.03008, 1.98907, 1.9179, 1.86916, 1.80686, 
      1.89057, 1.8693, 1.8615, 1.81276, 1.25726, 1.72066, 1.75852, 
      1.68129, 1.01656, 2.04999, 1.77068, -4.73643, -4.79763, -4.78251, 
      -4.76003, -4.9176, -4.74279, -4.98728, -4.96963, -4.21578, -4.21073, 
      -4.2797, -4.0338, -4.13588, -3.90354, -3.98899, -3.73241, -3.77705, 
      -3.66637, -3.57403, -3.54981, -3.39775, -3.23674, -3.21135, -3.13549, 
      -3.04044, -2.75852, -2.59817, -2.49531, -2.24252, -1.92726, -1.63032, 
      -0.95827, -1.28946, -2.3648, -2.37093, -1.966, -1.84397, -2.49355, 
      -2.52775, -5.31508, -5.17455, -5.06701, -5.0011, -4.91619, -4.80564, 
      -4.76309, -4.69903, -4.63178, -4.56894, -4.52929, -4.43313, -4.33866, 
      -4.28378, -4.21924, -4.18088, -4.12436, -4.05512, -3.99356, -3.9108, 
      -3.83678, -3.76069, -3.68691, -3.61607, -3.55028, -3.46728, -3.37619, 
      -3.28577, -3.15898, -3.04281, -2.92248, -2.79052, -2.64375, -2.47478, 
      -2.28517, -2.04776, -1.74824, -1.31333, -0.98662
    )
    
    dnorm_mod$obj$env$last.par.best <- c(
      delta_i = 0, logtau_i = log(0.5), logsigma_i = log(0.5), 
      epsilon_ti = epsilon_ti, nu_ti = nu_ti, logF_ti = logF_ti
    )
    
    dgmrf_mod$obj$env$last.par.best <- c(
      delta_i = 0, epsilon_ti = epsilon_ti, nu_ti = nu_ti, 
      beta = rep(0.5, 2), logF_ti = logF_ti
    )
    
    expect_equal(dnorm_mod$obj$fn(dnorm_mod$opt$par), dgmrf_mod$obj$fn(dgmrf_mod$opt$par))
    
  }
)


test_that(
  "fixing SEM parameters works", {
    
    # Some fixed, some not
    sem <- "
      eps_Pollock <-> eps_Pollock, 0, NA, 1, 
      eps_Arrowtooth <-> eps_Arrowtooth, 0, tau_atf, 0.1
    "
    
    dgmrf_mod <- ecostate(
      taxa = taxa, years = years, catch = ebs$Catch, biomass = ebs$Survey, 
      PB = ebs$P_over_B, QB = ebs$Q_over_B, B = B, DC = ebs$Diet_proportions, 
      EE = EE, X = X, type = type, U = U, fit_B0 = "Pollock",
      control = ecostate_control(nlminb_loops = 0, getsd = FALSE), 
      sem = sem
    )
    
    expect_equal(dgmrf_mod$sem$Type, c("Fixed", "Estimated"))
    expect_equal(dgmrf_mod$sem$Estimate, c(1, 0.1))
    expect_equal(dgmrf_mod$tmb_inputs$p$beta, c(tau_atf = 0.1))
    expect_equal(dgmrf_mod$tmb_inputs$map$beta, factor(1))
    
    # All fixed
    sem <- "
      eps_Pollock <-> eps_Pollock, 0, NA, 1, 
      eps_Arrowtooth <-> eps_Arrowtooth, 0, NA, 0.1
    "
    
    dgmrf_mod <- ecostate(
      taxa = taxa, years = years, catch = ebs$Catch, biomass = ebs$Survey, 
      PB = ebs$P_over_B, QB = ebs$Q_over_B, B = B, DC = ebs$Diet_proportions, 
      EE = EE, X = X, type = type, U = U, fit_B0 = "Pollock",
      control = ecostate_control(nlminb_loops = 0, getsd = FALSE), 
      sem = sem
    )
    
    expect_equal(dgmrf_mod$sem$Type, rep("Fixed", 2))
    expect_equal(dgmrf_mod$sem$Estimate, c(1, 0.1))
    expect_equal(dgmrf_mod$tmb_inputs$p$beta, numeric(0))
    expect_equal(dgmrf_mod$tmb_inputs$map$beta, NULL)
    
    # Equality constraint
    sem <- "
      eps_Pollock <-> eps_Pollock, 0, NA, 1, 
      eps_Arrowtooth <-> eps_Arrowtooth, 0, tau, 0.1
      eps_Cod <-> eps_Cod, 0, tau, 0.1
    "
    
    dgmrf_mod <- ecostate(
      taxa = taxa, years = years, catch = ebs$Catch, biomass = ebs$Survey, 
      PB = ebs$P_over_B, QB = ebs$Q_over_B, B = B, DC = ebs$Diet_proportions, 
      EE = EE, X = X, type = type, U = U, fit_B0 = "Pollock",
      control = ecostate_control(nlminb_loops = 0, getsd = FALSE), 
      sem = sem
    )
    
    expect_equal(dgmrf_mod$sem$Type, c("Fixed" , "Estimated", "Estimated"))
    expect_equal(dgmrf_mod$sem$name, c(NA , "tau", "tau"))
    expect_equal(dgmrf_mod$sem$Estimate, c(1, 0.1, 0.1))
    expect_equal(dgmrf_mod$tmb_inputs$p$beta, c(tau = 0.1))
    expect_equal(dgmrf_mod$tmb_inputs$map$beta, factor(1))
    
  }
)

test_that(
  "SEM works with missing covariates", {
    
    sem <- "
      cold_pool -> eps_Pollock, 1, beta_cp_plk, 1,
      eps_Pollock <-> eps_Pollock, 0, tau_plk, 1
    "
    
    expect_error(ecostate(
      taxa = taxa, years = years, catch = ebs$Catch, biomass = ebs$Survey, 
      PB = ebs$P_over_B, QB = ebs$Q_over_B, B = B, DC = ebs$Diet_proportions, 
      EE = EE, X = X, type = type, U = U, fit_B0 = "Pollock",
      control = ecostate_control(nlminb_loops = 0, getsd = FALSE), 
      sem = sem
    ), "Some variable(s) in `sem` are not in `covariates`", fixed = TRUE) 
    
    dgmrf_mod <- ecostate(
      taxa = taxa, years = years, catch = ebs$Catch, biomass = ebs$Survey, 
      PB = ebs$P_over_B, QB = ebs$Q_over_B, B = B, DC = ebs$Diet_proportions, 
      EE = EE, X = X, type = type, U = U, fit_B0 = "Pollock",
      control = ecostate_control(nlminb_loops = 0, getsd = FALSE), 
      sem = sem, covariates = covariates
    )
    
    expect_equal(dgmrf_mod$tmb_inputs$p$covariates, covariates)
    expect_equal(dgmrf_mod$tmb_inputs$map$covariates, factor(rep(NA, nrow(covariates))))
    
    covariates[5:10,] <- NA
    
    dgmrf_mod <- ecostate(
      taxa = taxa, years = years, catch = ebs$Catch, biomass = ebs$Survey, 
      PB = ebs$P_over_B, QB = ebs$Q_over_B, B = B, DC = ebs$Diet_proportions, 
      EE = EE, X = X, type = type, U = U, fit_B0 = "Pollock",
      control = ecostate_control(nlminb_loops = 0, getsd = FALSE), 
      sem = sem, covariates = covariates
    )
    
    map_expected <- seq_len(nrow(covariates))
    map_expected <- factor(ifelse(is.na(covariates), map_expected, NA))
    expect_equal(dgmrf_mod$tmb_inputs$p$covariates, ifelse(is.na(covariates), 0, covariates))
    expect_equal(dgmrf_mod$tmb_inputs$map$covariates, map_expected)
    
  }
)

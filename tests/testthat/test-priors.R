
test_that(
  "Parsing density list works", {
    
    p <- list(
      x = setNames((-2):2, letters[1:5]), 
      y = 0,
      z = matrix(1:9, ncol = 3, nrow = 3)
    )
    
    # Global evaluation
    lprior <- sum(dnorm(p$x[3:5], log = TRUE, mean = 1, sd = 2)) + 
      dexp(exp(p$y), log = TRUE) + 
      sum(dlnorm(diag(p$z), meanlog = 1, log = TRUE))
    
    # List implementation
    prior_list <- list(
      x[c("c", "d", "e")] ~ dnorm(mean = 1, sd = 2), 
      exp(y) ~ dexp(), 
      diag(z) ~ dlnorm(meanlog = 1)
    )
    
    # Function implementation
    lprior_fn <- function(p) {
      lprior <- sum(dnorm(p$x[c("c", "d", "e")], mean = 1, sd = 2, log = TRUE))
      lprior <- lprior + dexp(exp(p$y), log = TRUE)
      lprior <- lprior + sum(dlnorm(diag(p$z), meanlog = 1, log = TRUE))
      lprior
    }
    
    expect_equal(evaluate_prior(prior_list, p), lprior)
    expect_equal(evaluate_prior(lprior_fn, p), lprior)
    expect_error(evaluate_prior(append(prior_list, x["f"] ~ dnorm()), p), "Problem with prior x")
    expect_error(evaluate_prior(append(prior_list, w ~ dnorm(sd = 2)), p), "Problem with prior w")
    
  }
)

test_that(
  "Parsing SEM priors works", {
    
    p <- list(
      beta = c(tau_plk = 0.1, sigma_atf = 0.1),
      DC_ij = structure(
        c(0.128, 0.001, 0.001, 0.872, 0.004, 0, 0.902, 0.002, 0.018), dim = c(3L, 3L), 
        dimnames = list(Prey = c("Pollock", "Arrowtooth", "Cod"), 
                        Predator = c("Pollock", "Arrowtooth", "Cod"))
      ), 
      epsilon_ti = matrix(0, nrow = 10, ncol = 3), 
      nu_ti =  matrix(0, nrow = 10, ncol = 3)
    )
    colnames(p$epsilon_ti) <- colnames(p$nu_ti) <- c("Pollock", "Arrowtooth", "Cod")
    
    lprior <- dnorm(log(abs(p$beta["tau_plk"])), log = TRUE) + 
      dlnorm(abs(p$beta["sigma_atf"]), sdlog = 0.1, log = TRUE) + 
      dnorm(p$epsilon_ti[1, "Pollock"], log = TRUE) + 
      sum(dnorm(qlogis(p$DC_ij[,"Pollock"]), mean = 0, sd = 1, log = TRUE))
    
    names(lprior) <- NULL
    
    prior_list <- list(
      log(abs(tau_plk)) ~ dnorm(), 
      abs(sigma_atf) ~ dlnorm(sdlog = 0.1), 
      epsilon_ti[1, "Pollock"] ~ dnorm(),
      qlogis(DC_ij[, "Pollock"]) ~ dnorm()
    )
    
    expect_equal(evaluate_prior(prior_list, p), lprior)
    
    p$beta["DC_ij"] <- 2
    expect_error(evaluate_prior(prior_list, p), "Parameter name DC_ij is reserved")
    
  }
)
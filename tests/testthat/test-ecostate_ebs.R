

test_that(
  "Published EBS biomass-dynamics model works", {

  data( eastern_bering_sea )
  attach( eastern_bering_sea )

  # Run model
  out = ecostate( taxa = taxa,
                  years = years,
                  catch = Catch,
                  biomass = Survey,
                  PB = P_over_B,
                  QB = Q_over_B,
                  DC = Diet_proportions,
                  B = fittedB_i,
                  EE = EE_i,
                  U = U_i,
                  type = type_i,
                  X = X_ij,
                  fit_B = fit_B,
                  fit_Q = fit_Q,
                  fit_eps = fit_eps,
                  fit_B0 = fit_B0,
                  fit_PB = fit_PB,
                  fit_QB = fit_QB,
                  log_prior = log_prior,
                  control = ecostate_control( n_steps = n_steps,
                                              getsd = FALSE,
                                              nlminb_loops = 0,
                                              silent = TRUE,
                                              start_tau = 0.01,
                                              tmbad.sparse_hessian_compress = 1,
                                              trace = 1 ) )

  #
  newval = out$obj$fn(opt_MLE$par)
  expect_equal( as.numeric(newval), opt_MLE$obj, tolerance = 0.01 )

  # Re-run from starting value of fixed and random effects given MLE
  tmb_par = out$obj$env$parList()
  out = ecostate( taxa = taxa,
                  years = years,
                  catch = Catch,
                  biomass = Survey,
                  PB = P_over_B,
                  QB = Q_over_B,
                  DC = Diet_proportions,
                  B = fittedB_i,
                  EE = EE_i,
                  U = U_i,
                  type = type_i,
                  X = X_ij,
                  fit_B = fit_B,
                  fit_Q = fit_Q,
                  fit_eps = fit_eps,
                  fit_B0 = fit_B0,
                  fit_PB = fit_PB,
                  fit_QB = fit_QB,
                  log_prior = log_prior,
                  control = ecostate_control( n_steps = n_steps,
                                              getsd = FALSE,
                                              nlminb_loops = 1,
                                              tmb_par = tmb_par,
                                              silent = TRUE,
                                              start_tau = 0.01,
                                              tmbad.sparse_hessian_compress = 1,
                                              trace = 1 ) )

  # Check methods
  print(out)
  logLik(out)

  # Check simulators
  set.seed(101)
  sim1 = out$simulator( simulate_random = FALSE )
  set.seed(101)
  sim2 = out$simulator( simulate_random = FALSE )
  set.seed(101)
  sim3 = out$simulator( simulate_random = TRUE )
  set.seed(101)
  sim4 = out$simulator( simulate_random = TRUE )
  expect_true( identical(sim1, sim2) )
  expect_true( identical(sim3, sim4) )
  expect_false( identical(sim1, sim3) )

  # Plots
  library(ggplot2)
  ggplot_foodweb = plot_foodweb( out$rep$out_initial$Qe_ij,
                xtracer_i = ifelse(taxa=="Chloro",1,0),
                B_i = out$rep$out_initial$B_i,
                type_i = type_i )

})




#' @title Plot foodweb
#'
#' @description 
#' Plot consumption as a directed graph including all taxa (vertices) 
#' and biomass consumed (arrows).  Taxa are located using tracers,
#' where by default the y-axis is trophic level.  #' 
#'
#' @inheritParams compute_tracer
#'
#' @param Q_ij Consumption of each prey i by predator j in units biomass.
#' @param xtracer_i tracer to use when computing x-axis values
#' @param ytracer_i tracer to use when computing y-axis values
#' @param B_i biomass to use when weighting taxa in plot
#' @param taxa_labels character vector of labels to use for each taxon
#' @param xloc x-axis location (overrides calculation using \code{xtracer_i})
#' @param yloc y-axis location (overrides calculation using \code{ytracer_i})
#' @param rescale whether to rescale flow and biomass as ratio relative to lowest value
#'        such that values are relative rather than absolute units
#' @param ... additional arguments to \code{ggnetwork::ggnetwork()}
#'
#' @importFrom ggplot2 .data
#' 
#' @details
#' Trophic level \eqn{l_i} for each predator \eqn{i} is defined as:
#'
#' \deqn{ \mathbf{l - 1 = l Q^*} }
#'
#' where \eqn{\mathbf{Q*}} is the proportion consumption for each predator (column)
#' of different prey (rows).  We identify primary producers as any taxa with no
#' consumption (a column of 0s), and assign them as the first trophic level.
#'
#' @return
#' invisibly return \code{ggplot} object for foodweb
#'
#' @examples
#' \dontrun{
#' # Load data
#' data(eastern_bering_sea)
#' 
#' # Simple model estimating only Pollock biomass
#' model <- with(eastern_bering_sea, ecostate(
#'   taxa = taxa, 
#'   years = years, 
#'   catch = Catch, 
#'   biomass = Survey, 
#'   type = type_i, 
#'   fit_B = "Pollock",
#'   B = B_i, 
#'   PB = P_over_B, 
#'   QB = Q_over_B,
#'   EE = EE_i, 
#'   DC = Diet_proportions,
#'   U = U_i, 
#'   X = X_ij
#' ))
#' 
#' # Plot using consumption at equilibrium Qe_ij
#' # X-axis is proportion of production from Krill (pelagic production)
#' plot_foodweb(
#'   model$rep$out_initial$Qe_ij, 
#'   xtracer_i = ifelse(eastern_bering_sea$taxa == "Krill", 1, 0), 
#'   B_i = model$rep$out_initial$B_i, 
#'   type_i = eastern_bering_sea$type_i,
#'   taxa_labels = eastern_bering_sea$taxa
#' )
#' }
#'
#' @export
plot_foodweb <-
function( Q_ij,
          type_i,
          xtracer_i,
          ytracer_i = rep(1,nrow(Q_ij)),
          B_i = rep( 1, nrow(Q_ij)), 
          taxa_labels = letters[seq_len(nrow(Q_ij))],
          xloc,
          yloc,
          rescale = TRUE,
          ... ){

  #
  if(missing(yloc)){
    yloc = compute_tracer( Q_ij, 
                           inverse_method = "Standard", 
                           tracer_i = ytracer_i,
                           type_i = type_i )[1,]
  }
  if(missing(xloc)){
    xloc = compute_tracer( Q_ij, 
                           inverse_method = "Standard", 
                           tracer_i = xtracer_i,
                           type_i = type_i )[1,]
  }

                
  #
  layout = cbind( x = xloc, 
                  y = yloc )
  dimnames(Q_ij) = list( taxa_labels, taxa_labels )
  #rownames(layout) = taxa_labels                       
  graph = igraph::graph_from_adjacency_matrix( Q_ij, 
                                       weighted = TRUE )
  #plot( graph,
  #           edge.width = log(E(graph)$weight),
  #           vertex.size = log(B_i/min(B_i)) + 10,
  #           type = "width", 
  #           text_size = 4,
  #           #arrow = grid::arrow(type='closed', 18, grid::unit(10,'points')),
  #           layout = layout,
  #           show.legend = FALSE, 
  #           curvature = 0 )
  #axis(1)
  #axis(2)
  #return( list('graph' = graph,
  #             'layout' = layout))

  g = ggnetwork::ggnetwork( x = graph, 
                 layout = layout,
                 scale = FALSE, 
                ...)
  #g = ggnetwork::ggnetwork( Q_ij, 
  #                          layout = layout,
  #                          weighted = TRUE )
  #g$log_flow = log(g$weight / min(g$weight,na.rm=TRUE))
  g$flow = g$weight
  if(rescale==TRUE) g$flow = g$flow / min(g$flow,na.rm=TRUE)
  g$mass = rep(NA, nrow(g))

  # Only fill in mass for rows without a flow, which represent node masses
  #g$log_mass[which(is.na(g$log_flow))] = log(B_i / min(B_i,na.rm=TRUE))
  if( length(which(is.na(g$flow))) != length(B_i) ) stop("Check plotting code")
  g$mass[which(is.na(g$flow))] = B_i
  if(rescale==TRUE) g$mass = g$mass / min(g$mass,na.rm=TRUE)

  p = ggplot2::ggplot(g, ggplot2::aes(x=.data$x, y=.data$y, xend=.data$xend, yend=.data$yend) ) +
    ggnetwork::geom_edges( ggplot2::aes(colour=log(.data$flow)) ) +  #
    ggnetwork::geom_nodes( ggplot2::aes(size=log(.data$mass) ) ) +  #
    ggnetwork::geom_nodetext( ggplot2::aes(label=.data$name), fontface="bold", col="red")  #
  print(p)
  return(invisible(p))
}

#' @title Plot biomass time series from fitted ecostate model
#'
#' @param model Fitted model object returned by `ecostate`
#' @param taxa Taxon names to plot biomass for. If missing, defaults to all modeled taxa.
#' @param observed Plot survey biomass for available taxa? Defaults to TRUE.
#' @param interval Coverage of plotted confidence bands. Can be a numeric value between 0 and 1, 
#'   `NA` (or `FALSE`) for no interval, or `"SE"` for a +/- 1 SE interval. Defaults to 0.95.
#' @param q_adj Where applicable, should estimates adjusted for catchability be plotted? Useful for 
#'   assessing fit to survey data.
#' @param dist Distribution assumption for biomass uncertainty. Either `"lognormal"` (default) or `"normal"`.
#'   
#' @importFrom stats reshape qlnorm qnorm
#' @importFrom ggplot2 .data
#'
#' @return invisibly return \code{ggplot} object for biomass time series
#' @export
plot_timeseries <- function(model, taxa, observed = TRUE, interval = 0.95, q_adj = TRUE, dist = c("lognormal", "normal")) {
  
  dist <- match.arg(dist)
  years <- union(model$internal$years, model$internal$extra_years)
  
  if (missing(taxa)) {
    taxa <- model$internal$taxa
  } else if (!all(taxa %in% model$internal$taxa)) {
    stop(paste("taxa not in model:", paste(taxa[!(taxa %in% model$internal$taxa)], collapse = ", ")))
  }
  
  # Estimated biomass
  if (is.null(model$derived$Est$B_ti)) stop("No derived biomass estimates to plot")
  B_ti <- list(est = model$derived$Est$B_ti)
  
  # Determine if we should plot confidence interval and check type
  plot_interval <- FALSE
  use_se <- FALSE
  
  if (!is.null(interval) && !(length(interval) == 1 && (is.na(interval) || isFALSE(interval)))) {
    if (identical(interval, "SE")) {
      plot_interval <- TRUE
      use_se <- TRUE
    } else if (is.numeric(interval) && length(interval) == 1) {
      if (!(interval > 0 && interval < 1)) stop("interval must be between 0 and 1")
      plot_interval <- TRUE
    } else {
      stop("interval must be a numeric value between 0 and 1, NA, or 'SE'")
    }
  }
  
  # Standard error
  if (plot_interval) {
    if (!all(is.na(model$derived$SE$B_ti))) {
      B_ti[["se"]] <- model$derived$SE$B_ti
    } else {
      message("`interval` is ignored when standard errors of derived quantities are missing")
      plot_interval <- FALSE
    }
  }
  
  # Observed (survey) biomass
  if (isTRUE(observed)) B_ti[["observed"]] <- model$internal$Bobs_ti
  
  # Estimated biomass, adjusted using catchability
  if (isTRUE(q_adj)) {
    qmat <- outer(rep(1, length(years)), exp(model$internal$parhat$logq_i))
    qmat[qmat == 1] <- NA
    B_ti[["est_q"]] = model$derived$Est$B_ti * qmat
  }
  
  # Convert all to long format
  B_ti <- B_ti |> 
    lapply(as.data.frame) |> 
    lapply(function(x) {colnames(x) <- model$internal$taxa; x}) |> 
    lapply(function(x) {x$year <- years; x}) |>  
    lapply(reshape, direction = "long", varying = model$internal$taxa,
      v.names = "biomass", times = model$internal$taxa, timevar = "Taxon"
    ) 
  
  # Join data frames
  for (i in seq_along(B_ti)) names(B_ti[[i]])[names(B_ti[[i]]) == "biomass"] <- names(B_ti)[i]
  B_ti <- Reduce(function(x, y) merge(x, y, c("year", "Taxon", "id")), B_ti)
  
  # Subset to specified taxa
  B_ti <- B_ti[B_ti$Taxon %in% taxa, ]
  
  # Add confidence interval
  if (plot_interval) {
    
    if (use_se) {
      if (dist == "lognormal") {
        B_ti$se_log <- sqrt(log(1 + ((B_ti$se^2) / (B_ti$est^2))))
        B_ti$est_log <- log(B_ti$est)
        B_ti$lower <- exp(B_ti$est_log - B_ti$se_log)
        B_ti$upper <- exp(B_ti$est_log + B_ti$se_log)
      } else {
        B_ti$lower <- B_ti$est - B_ti$se
        B_ti$upper <- B_ti$est + B_ti$se
      }
    } else {
      if (dist == "lognormal") {
        B_ti$se_log <- sqrt(log(1 + ((B_ti$se^2) / (B_ti$est^2))))
        B_ti$est_log <- log(B_ti$est)
        B_ti$lower <- qlnorm((1 - interval) / 2, meanlog = B_ti$est_log, sdlog = B_ti$se_log)
        B_ti$upper <- qlnorm(1 - (1 - interval) / 2, meanlog = B_ti$est_log, sdlog = B_ti$se_log)
      } else {
        z_val <- qnorm(1 - (1 - interval) / 2)
        B_ti$lower <- B_ti$est - z_val * B_ti$se
        B_ti$upper <- B_ti$est + z_val * B_ti$se
      }
    }
    
  }
  
  colors <- c("Bhat" = "black")
  
  # Plot estimate and interval
  p <- B_ti |> 
    ggplot2::ggplot(ggplot2::aes(.data$year)) + 
    {if (plot_interval) ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), alpha = 0.25, na.rm = TRUE)} + 
    ggplot2::geom_line(ggplot2::aes(y = .data$est, color = "Bhat"), na.rm = TRUE) + 
    ggplot2::facet_wrap(~.data$Taxon, scales = "free_y") + 
    ggplot2::labs(y = "biomass", color = "source")
  
  # Plot estimate adjusted for catchability
  if (isTRUE(q_adj)) {
    colors <- c(colors, "Bhat * q" = "blue")
    p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$est_q, color = "Bhat * q"), na.rm = TRUE) 
  }
  
  # Plot survey data
  if (isTRUE(observed)) {
    colors <- c(colors, "Survey" = "black")
    p <- p + ggplot2::geom_point(ggplot2::aes(y = .data$observed, color = "Survey"), na.rm = TRUE)
  }
  
  # Set colors
  p <- p + ggplot2::scale_color_manual(values = colors)
  
  print(p)
  return(invisible(p))
  
}


#' @title Plot fits to age composition data
#'
#' @param model Fitted model returned by `ecostate`
#' @param stgroup Stanza group to plot age composition fits for. Must be length 1.
#'        Defaults to the first stanza group in \code{settings$unique_stanza_groups}.
#' @param years Years to plot. Defaults to all years with age composition data
#' @param track_cohorts Whether to color observed agecomp data by cohort. Defaults to TRUE
#' @param cutoff Minimum proportion for plotting an age bin. Specifying a cutoff may 
#'        reduce the number of cohorts and make cohort colors easier to discern.
#'
#' @importFrom stats reshape
#' @importFrom ggplot2 ggplot aes .data
#'
#' @returns A ggplot object
#' @export
plot_agecomp <- function(model, stgroup, years, track_cohorts = TRUE, cutoff = 0) {
  
  if (missing(stgroup)) stgroup <- model$internal$settings$unique_stanza_groups[1]
  
  obs <- model$internal$Nobs_ta_g2[[stgroup]]
  if (is.null(obs)) stop("No age composition data available for this stanza group")
  if (missing(years)) years <- as.integer(rownames(obs))
  
  # Convert to proportions
  est <- t(apply(model$rep$Nexp_ta_g2[[stgroup]], 1, function(x) x/sum(x, na.rm = TRUE)))
  est <- as.data.frame(est)
  est$year <- as.integer(rownames(est))
  
  est <- reshape(
    est, direction = "long", 
    varying = colnames(est)[colnames(est) != "year"], 
    v.names = "estimated", timevar = "age", idvar = "year"
  )
  
  # Convert to proportions
  obs <- t(apply(obs, 1, function(x) x/sum(x, na.rm = TRUE)))
  
  obs <- as.data.frame(obs)
  obs$year <- as.integer(rownames(obs))
  
  obs <- reshape(
    obs, direction = "long",
    varying = colnames(obs)[colnames(obs) != "year"], 
    v.names = "observed", timevar = "age", idvar = "year"
  )
  
  acomp <- merge(est, obs, by = c("year", "age"))
  
  # Track cohorts
  if (isTRUE(track_cohorts)) {
    
    cohorts <- outer(
      seq_len(nrow(model$rep$Nexp_ta_g2[[stgroup]])), 
      seq_len(ncol(model$rep$Nexp_ta_g2[[stgroup]])), 
      FUN = "-"
    )
    
    cohorts[] <- as.integer(as.factor(cohorts[]))
    
    dimnames(cohorts) <- dimnames(model$rep$Nexp_ta_g2[[stgroup]])
    
    cohorts <- as.data.frame(cohorts)
    cohorts$year <- as.integer(rownames(cohorts))
    
    cohorts <- reshape(
      cohorts, direction = "long",
      varying = colnames(cohorts)[colnames(cohorts) != "year"], 
      v.names = "cohort", timevar = "age", idvar = "year"
    )
    
    acomp <- merge(acomp, cohorts, by = c("year", "age"))
    
  }
  
  acomp <- acomp[acomp$year %in% years,]
  
  acomp$observed[acomp$observed < cutoff] <- NA
  
  p <- acomp |> ggplot2::ggplot(ggplot2::aes(.data$age, .data$observed))
  
  if (isTRUE(track_cohorts)) {
    
    # Copy of viridis::turbo palette to reduce dependencies
    pal <- grDevices::colorRampPalette(c(
      "#30123BFF", "#4454C4FF", "#4490FEFF", "#1FC8DEFF", "#29EFA2FF", "#7DFF56FF", 
      "#C1F334FF", "#F1CA3AFF", "#FE922AFF", "#EA4F0DFF", "#BE2102FF", "#7A0403FF"
     ))
    
    colors <- pal(length(unique(acomp$cohort)))
    
    colors_dark <- vapply(
      colors, \(x) grDevices::colorRampPalette(c(x, "black"))(5)[3], 
      character(1), USE.NAMES = FALSE
    )
    
    p <- p + 
      ggplot2::geom_bar(
        ggplot2::aes(fill = factor(.data$cohort), color = factor(.data$cohort)), 
        stat = "identity", show.legend = FALSE, na.rm = TRUE, width = 1
      ) + 
      ggplot2::scale_fill_manual(values = colors) + 
      ggplot2::scale_color_manual(values = colors_dark)
    
  } else {
    
    p <- p + ggplot2::geom_bar(
      stat = "identity", show.legend = FALSE, width = 1, 
      na.rm = TRUE, fill = "grey40", color = "grey20"
    )
    
  }
  
  p <- p + 
    ggplot2::geom_line(ggplot2::aes(y = .data$estimated)) +
    ggplot2::geom_point(ggplot2::aes(y = .data$estimated)) + 
    ggplot2::facet_wrap(~.data$year, dir = "v") + 
    ggplot2::ylab("proportion")
  
  print(p)
  return(invisible(p))
  
}

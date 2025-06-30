

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
          rescale = TRUE ){

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
                 scale = FALSE )
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

  # https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
  x = y = xend = yend = name = flow = mass = NULL
  p = ggplot2::ggplot(g, ggplot2::aes(x=x, y=y, xend=xend, yend=yend) ) +
    ggnetwork::geom_edges( ggplot2::aes(colour=log(flow)) ) +  #
    ggnetwork::geom_nodes( ggplot2::aes(size=log(mass) ) ) +  #
    ggnetwork::geom_nodetext( ggplot2::aes(label=name), fontface="bold", col="red")  #
  print(p)
  return(invisible(p))
}

#' Plot biomass time series from fitted ecostate model
#'
#' @param fit Fitted model object returned by `ecostate`
#' @param taxa Taxon names to plot biomass for. If missing, defaults to all modeled taxa.
#' @param observed Plot survey biomass for available taxa? Defaults to TRUE.
#' @param interval Coverage of plotted confidence bands. Defaults to 0.95 (95%), omitted if FALSE.
#' @param q_adj Where applicable, should estimates adjusted for catchability be plotted? Useful for 
#'   assessing fit to survey data.
#'
#' @return invisibly return \code{ggplot} object for biomass time series
#' @export
plot_timeseries <- function(fit, taxa, observed = TRUE, interval = 0.95, q_adj = TRUE) {
  
  if (missing(taxa)) {
    taxa <- fit$internal$taxa
  } else if (!all(taxa %in% fit$internal$taxa)) {
    stop(paste("taxa not in model:", paste(taxa[!(taxa %in% fit$internal$taxa)], collapse = ", ")))
  }
  
  # Estimated biomass
  if (is.null(fit$derived$Est$B_ti)) stop("No derived biomass estimates to plot")
  B_ti <- list(est = fit$derived$Est$B_ti)
  
  # Standard error
  if (!isFALSE(interval)) {
    if (!all(is.na(fit$derived$SE$B_ti))) {
      B_ti[["se"]] <- fit$derived$SE$B_ti
    } else {
      message("`interval` is ignored when standard errors of derived quantities are missing")
      interval <- FALSE
    }
  }
  
  # Observed (survey) biomass
  if (isTRUE(observed)) B_ti[["observed"]] <- fit$internal$Bobs_ti
  
  # Estimated biomass, adjusted using catchability
  if (isTRUE(q_adj)) {
    qmat <- outer(rep(1, length(fit$internal$years)), exp(fit$internal$parhat$logq_i))
    qmat[qmat == 1] <- NA
    B_ti[["est_q"]] = fit$derived$Est$B_ti * qmat
  }
  
  # Convert all to long format
  B_ti <- B_ti |> 
    lapply(as.data.frame) |> 
    lapply(function(x) {colnames(x) <- fit$internal$taxa; x}) |> 
    lapply(function(x) {x$year <- fit$internal$years; x}) |>  
    lapply(reshape, direction = "long", varying = fit$internal$taxa,
      v.names = "biomass", times = fit$internal$taxa, timevar = "Taxon"
    ) 
  
  # Join data frames
  for (i in seq_along(B_ti)) names(B_ti[[i]])[names(B_ti[[i]]) == "biomass"] <- names(B_ti)[i]
  B_ti <- Reduce(function(x, y) merge(x, y, c("year", "Taxon", "id")), B_ti)
  
  # Subset to specified taxa
  B_ti <- B_ti[B_ti$Taxon %in% taxa, ]
  
  # Add confidence interval
  if (!isFALSE(interval)) {
    
    checkmate::assert_number(interval)
    if (!(interval > 0 & interval < 1)) stop("interval must be between 0 and 1")
    
    B_ti$se_log <- sqrt(log(1 + ((B_ti$se^2) / (B_ti$est^2))))
    B_ti$est_log <- log(B_ti$est) - (0.5 * B_ti$se_log^2)
    B_ti$lower <- qlnorm((1 - interval) / 2, meanlog = B_ti$est_log, sdlog = B_ti$se_log)
    B_ti$upper <- qlnorm(1 - (1 - interval) / 2, meanlog = B_ti$est_log, sdlog = B_ti$se_log)
    
  }
  
  colors <- c("Bhat" = "black")
  
  # Plot estimate and interval
  p <- B_ti |> 
    ggplot2::ggplot(ggplot2::aes(year)) + 
    {if (!isFALSE(interval)) ggplot2::geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, na.rm = TRUE)} + 
    ggplot2::geom_line(ggplot2::aes(y = est, color = "Bhat"), na.rm = TRUE) + 
    ggplot2::facet_wrap(~Taxon, scales = "free_y") + 
    ggplot2::labs(y = "biomass", color = "source")
  
  # Plot estimate adjusted for catchability
  if (isTRUE(q_adj)) {
    colors <- c(colors, "Bhat * q" = "blue")
    p <- p + ggplot2::geom_line(ggplot2::aes(y = est_q, color = "Bhat * q"), na.rm = TRUE) 
  }
  
  # Plot survey data
  if (isTRUE(observed)) {
    colors <- c(colors, "Survey" = "black")
    p <- p + ggplot2::geom_point(ggplot2::aes(y = observed, color = "Survey"), na.rm = TRUE)
  }
  
  # Set colors
  p <- p + ggplot2::scale_color_manual(values = colors)
  
  print(p)
  return(invisible(p))
  
}

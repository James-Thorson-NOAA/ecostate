% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_equilibrium.R
\name{add_equilibrium}
\alias{add_equilibrium}
\title{Compute equilibrium values}
\usage{
add_equilibrium(ecoparams, scale_solver, noB_i, type_i)
}
\arguments{
\item{ecoparams}{list of parameters}

\item{scale_solver}{Whether to solve for ecotrophic efficiency EE given biomass B
(\code{scale_solver="simple"}) or solve for a combination of EE and B values}

\item{noB_i}{Boolean vector indicating which taxa have no B value}

\item{type_i}{character vector indicating whether a taxon is "hetero", "auto", or
"detritus"}
}
\value{
the list of parameters with missing values in \code{ecoparams$B_i} and/or
\code{ecoparams$EE_i} filled in, as well as additional values \code{Qe_ij},
\code{m0_i}, and \code{GE_i}
}
\description{
Compute equilibrium values
}
\details{
Replaces NA values in ecotrophic efficiency and/or biomass with equilibrium
solution, and then calculates equilibrium consumption, natural mortality,
and other rates.
}

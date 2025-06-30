
#' @title 
#' eastern Bering Sea ecosystem data
#'
#' @description 
#' Data used to demonstrate a Model of Intermediate Complexity (MICE)
#' for the eastern Bering Sea.  
#' `data(eastern_bering_sea)` loads a list that includes four components:
#' * `Survey` is a long-form data-frame with three columns, providing the Year,
#'    Mass (in relative units for most taxa, and million metric tons for Pollock,
#'    Cod, Arrowtooth, and NFS), and Taxon for each year with available data
#' * `Catch` is a long-form data-frame with three columns, providing the Year,
#'    Mass (in million metric tons), and Taxon for each year with available data
#' * `P_over_B` is a numeric vector with the unitless ratio of biomass production to
#'   population biomass for each taxon (as vector names)
#' * `Q_over_B` is a numeric vector with the unitless ratio of biomass consumption to
#'   population biomass for each taxon (as vector names)
#' * `Diet_proportions` is a numeric matrix where each column lists the
#'   proportion of biomass consumed for each predator (column names)
#'   that is provided by each prey (row names)
#'
#' @details
#' The data compiled come from a variety of sources:
#' * Northern fur seal (NFS) survey is an absolute index, corrected for proportion
#'   of time spent in the eastern Bering Sea.  NFS QB is developed from a bioenergetic 
#'   model and also corrected for seasonal residency. Both are provided by 
#'   Elizabeth McHuron. It is post-processed in a variety of ways, and not
#'   to be treated as an index of abundance for NFS for other uses.
#' * Pollock, cod, and arrowtooth surveys are from a bottom trawl survey, and 
#'   cod and arrowtooth are treated as an absolute index.
#' * Copepod and other zooplankton are from an oblique tow bongo net survey, 
#'   with data provided by Dave Kimmel.  It is then post-processed to account
#'   for spatially and seaonally imbalanced data.
#' * Other P_over_B, Q_over_B and Diet_proportions values 
#'   are derived from Rpath models, provided by Andy Whitehouse.
#' * Primary producers is an annual index of relative biomass, developed from monthly
#'   satellite measurements and provided by Jens Nielsen.
#' See Thorson et al. (2025) for more details regarding data standardization
#'   and sources
#'
#' @name eastern_bering_sea
#' @docType data
#' @usage data(eastern_bering_sea)
#' @keywords data
NULL

#' @title
#' eastern Gulf of Alaska data
#'
#' @description
#' Data used to demonstrate size-structured modelling including bottom-up
#'   linkages.
#' `data(gulf_of_alaska)` loads a list that includes four components:
#' * `years` the set of years modeled
#' * `taxa` the set of taxa (including each stanza for age-structured populations)
#' * `type` is a character vector indicating for each taxon (as vector names)
#'    whether it is a heterotroph, autotroph, or detritus pool
#' * `biomass_data` is a long-form data-frame with three columns, providing the Year,
#'    Mass (in relative units for copepods and euphausiids, and million metric
#'    tons for Pollock and sablefish), and Taxon for each year with available data
#' * `catch_data` is a long-form data-frame with three columns, providing the Year,
#'    Mass (in million metric tons), and Taxon for each year with available data
#' * `P_over_B` is a numeric vector with the unitless ratio of biomass production to
#'   population biomass for each taxon (as vector names)
#' * `Q_over_B` is a numeric vector with the unitless ratio of biomass consumption to
#'   population biomass for each taxon (as vector names)
#' * `EE` is a numeric vector with the ecotrophic efficiency for each `taxa`, in
#'   this case provided only for the detritus pool and determining its turn-over rate
#' * `Diet_proportions` is a numeric matrix where each column lists the
#'   proportion of biomass consumed for each predator (column)
#'   that is provided by each prey (row)
#' * `agecomp_data` is a named list with age-composition samples for any
#'   age-structured population (named list element), where each list element
#'   is a matrix with named rows indicating years, and columns indicating integer
#'   ages
#' * `stanza_groups` is a character vector with names matching taxa, and elements
#'   indicating multi-stanza groups (i.e., age-structured populations)
#' * `Amax` is a numeric vector indicating the maximum age (in years) for each 
#'   stanza as indicated in `names(stanza_groups)`.
#' * `Leading` is a Boolean vector indicating which stanza is used for estimating
#'   biomass (and therefore recruitment) for each `names(stanza_groups)`.
#' * `K` is a numeric vector with the generalized von Bertalanffy growth coefficient
#'   (in weight not length) for each `unique(stanza_groups)`
#' * `d` is a numeric vector with the generalized von Bertalanffy growth coefficient
#'   (in weight not length) for each `unique(stanza_groups)`
#' * `Amat` is a numeric vector with the age at 50% maturity (in years)
#'   for each `unique(stanza_groups)`. If provided, it is used to back-calculate
#'   and replace the value for `Wmat`
#' * `Wmat` is a numeric vector with the weight at 50% maturity (relative to W_inf)
#'   for each `unique(stanza_groups)`.
#' * `Wmatslope` is a numeric vector with the slope of the logistic maturity at weight
#'   (relative to W_inf) for each `unique(stanza_groups)`
#'
#' @details
#' The data compiled come from a variety of sources:
#' * Sablefish biomass index and age-comps are from a cooperative longline survey, and
#'   the index is treated as a relative index.
#' * Pollock biomass index and age-comps are from a bottom trawl survey, and
#'   the index is treated as a relative index with a prior on the catchability coefficient.
#' * Copepod biomass is from an oblique tow bongo net survey,
#'   with data provided by Dave Kimmel.
#' * Euphausiid biomass is from the Seward Long-Term Ecological Research program,
#'   provided by Russ Hopcroft to the Ecological Status Report.
#' * Other P_over_B, Q_over_B and Diet_proportions values
#'   are derived from Rpath models, provided by Bia Dias.
#' * Age-structed inputs are determined based on life-history information and
#'   biological samples.
#' See Thorson et al. (In review) for more details regarding data standardization
#'   and sources
#'
#' @name gulf_of_alaska
#' @docType data
#' @usage data(gulf_of_alaska)
#' @keywords data
NULL

#' @title
#' Full rpath inputs for eastern Bering Sea
#'
#' @description
#' All Rpath inputs from Whitehouse et al. 2021
#'
#' @name whitehouse_2021
#' @docType data
#' @usage data(whitehouse_2021)
#' @keywords data
NULL

#' @title
#' Output from optimizer for GOA stanza example
#'
#' @description
#' Objective function and MLE for coefficients in GOA stanza example
#'
#' @name goa_mle
#' @docType data
#' @usage data(goa_mle)
#' @keywords data
NULL


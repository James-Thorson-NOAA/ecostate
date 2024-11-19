# EcoState

[![Documentation](https://img.shields.io/badge/documentation-ecostate-orange.svg?colorB=E91E63)](https://james-thorson-noaa.github.io/ecostate/)

Package _ecostate_ fits a state-space mass-balance model _EcoState_ intended for aquatic ecosystems, using mass-balance equations matching those from Ecopath and dynamical equations matching Ecosim.  Unlike the Ecopath-with-Ecosim (_EwE_) package, _ecostate_ fits both biological parameters (e.g., equilibrium biomass and predator-prey vulnerability) and measurement parameters (e.g., catchability coefficients) via fit to time-series data.  _ecostate_ also estimates additional process errors representing nonstationarity in growth efficiency, ecotrophic efficient, migration, or other unmodeled processes.  These process errors allow biomass patterns to closely match available data, so that resulting consumption (and associated productivity and mortality rates) can accurately be conditioned upon any residual patterns.     

## Installation

_ecostate_ can be installed from GitHub using:

``` r
library(remotes)
install_github( "James-Thorson-NOAA/ecostate" )
```

Or to access vignettes from your R session, please instead use:

``` r
remotes::install_github( "James-Thorson-NOAA/ecostate",
                          build_vignettes = TRUE )
browseVignettes("ecostate")
```

# More details 

For more background please read:

Thorson, J.  Kristensen, K., Aydin, K., Gaichas, S., Kimmel, D.G., McHuron, E.A., Nielsen, J.N., Townsend, H., Whitehouse, G.A. EcoState:  Extending Ecopath with Ecosim to estimate biological parameters and process errors using RTMB and time-series data.  Pre-print URL: https://doi.org/10.32942/X2QK81 

or the updated version:

Thorson, J.  Kristensen, K., Aydin, K., Gaichas, S., Kimmel, D.G., McHuron, E.A., Nielsen, J.N., Townsend, H., Whitehouse, G.A. (In press) The benefits of hierarchical ecosystem models: demonstration using a new state-space mass-balance model EcoState.  Fish and Fisheries.  

# Create new package environment
.pkgenv <- new.env(parent = emptyenv())

# For unloading a package
.onUnload <- function(libpath) {}

# What to do on attaching the package
.onAttach <- function(libname, pkgname) {
  # packageStartupMessage("############################")
  # packageStartupMessage("Loading ibis package ...")
  # packageStartupMessage("############################")
  # if(!"INLA" %in% utils::installed.packages()[,1] ){
  #   packageStartupMessage("Installing INLA...")
  #   utils::install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")
  # }
  # Don't show constant rgdal warning
  options("rgdal_show_exportToProj4_warnings" = "none")

  # Set some default ibis options
  options('ibis.nthread' = parallel::detectCores() - 1)
  options('ibis.runparallel' = FALSE)
  options('ibis.setupmessages' = TRUE)
  # Known seed
  options('ibis.seed' = 19372)
  # Known engines
  options('ibis.engines' = c('GDB-Model','BART-Model',
                             'INLABRU-Model','BREG-Model','GLMNET-Model',
                             'INLA-Model','STAN-Model','XGBOOST-Model'))
  # Names of priors
  options('ibis.priors' = c('INLAPrior', 'BARTPrior', 'GDBPrior','GLMNETPrior',
                            'XGBPrior', 'BREGPrior', 'STANPrior'))
  # Use the future package for any options. Default is FALSE
  options('ibis.futurestrategy' = "multisession")
  options(doFuture.foreach.export = ".export-and-automatic-with-warning")

  # Other dependencies not directly added in DESCRIPTION (to minimize potential issues)
  options('ibis.dependencies' = c(
    "pdp", "scales", "biscale", "modEvA", "dplyr", "geodist", "geosphere", "progress",
    "glmnet", "glmnetUtils", "xgboost","BoomSpikeSlab", "INLA", "inlabru",
    "dbarts", "mboost", "rstan", "cmdstanr"
  ))

  # Set default corrrelation coefficient threshold for evaluating correlated predictors
  options('ibis.corPred' = 0.7)

  # Set default settings for pseudo-absence sampling
  options('ibis.pseudoabsence' = pseudoabs_settings() )

  # Set S2 use for SF to false owing to the multiple bugs and errors with
  # 29/06 To be changed later eventually
  suppressMessages( sf::sf_use_s2(FALSE) )
  # Use sf instead of sp where possible
  # https://r-spatial.org/r/2023/05/15/evolution4.html
  options("sp_evolution_status" = 2)
}

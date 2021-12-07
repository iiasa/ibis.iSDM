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
  options('ibis.runparallel' = TRUE)
  options('ibis.setupmessages' = TRUE)
  options('ibis.engines' = c('GDB-Model','BART-Model',
                             'INLABRU-Model',
                             'INLA-Model','STAN-Model','XGBOOST-Model'))
  # Set S2 use for SF to false owing to the multiple bugs and errors with
  # 29/06 To be changed later eventually
  suppressMessages( invisible( sf::sf_use_s2(FALSE) ) )
}

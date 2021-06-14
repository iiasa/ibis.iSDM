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

  # Set some default ibis options
  options('ibis.nthread' = parallel::detectCores() - 1)
  options('ibis.runparallel' = FALSE)
}

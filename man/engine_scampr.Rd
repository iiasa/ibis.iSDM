% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/engine_scampr.R
\name{engine_scampr}
\alias{engine_scampr}
\title{Engine for process models using scampr}
\usage{
engine_scampr(x, type = "response", dens = "posterior", maxit = 500, ...)
}
\arguments{
\item{x}{\code{\link[=distribution]{distribution()}} (i.e. \code{\linkS4class{BiodiversityDistribution}}) object.}

\item{type}{The mode used for creating (posterior or prior) predictions. Either stting
\code{"link"} or \code{"response"} (Default: \code{"response"}).}

\item{dens}{A \code{\link{character}} on how predictions are made, either from the \code{"posterior"} (Default)
or \code{"prior"}.}

\item{maxit}{A \code{\link{numeric}} on the number of iterations for the optimizer (Default: \code{500}).}

\item{...}{Other parameters passed on.}
}
\value{
An \link{Engine}.
}
\description{
Similar to others, this engine enables the fitting and prediction of
log-Gaussian Cox process (LGCP) and Inhomogeneous Poisson process (IPP) processes.
It uses the \code{scampr} package, which uses maximum likelihood estimation
fitted via \code{TMB} (Template Model Builder).

It also support the addition of spatial latent effects which can be added via
Gaussian fields and approximated by 'FRK' (Fixed Rank Kriging) and are
integrated out using either variational or Laplace approximation.

The main use case for this engine is as an alternative to \code{\link[=engine_inlabru]{engine_inlabru()}} and
\code{\link[=engine_inla]{engine_inla()}} for fitting iSDMs, e.g. those combining both presence-only
and presence-absence point occurrence data.
}
\details{
This engine may only be used to predict for one or two datasets at most. It
supports only presence-only PPMs and presence/absence Binary GLMs, or 'IDM'
(for an integrated data model).
}
\note{
\itemize{
\item The package can currently be installed from github directly only \code{"ElliotDovers/scampr"}
\item Presence-absence models in SCAMPR currently only support cloglog link functions!
}
}
\examples{
\dontrun{
# Load background
background <- terra::rast(system.file('extdata/europegrid_50km.tif',
package='ibis.iSDM',mustWork = TRUE))

# Add GLM as an engine
x <- distribution(background) |> engine_scampr()
}
}
\references{
\itemize{
\item Dovers, E., Popovic, G. C., & Warton, D. I. (2024). A fast method for fitting integrated species distribution models. Methods in Ecology and Evolution, 15(1), 191-203.
\item Dovers, E., Stoklosa, D., and Warton D. I. (2024). Fitting log-Gaussian Cox processes using generalized additive model software. The American Statistician, 1-17.
}
}
\seealso{
Other engine: 
\code{\link{engine_bart}()},
\code{\link{engine_breg}()},
\code{\link{engine_gdb}()},
\code{\link{engine_glm}()},
\code{\link{engine_glmnet}()},
\code{\link{engine_inla}()},
\code{\link{engine_inlabru}()},
\code{\link{engine_stan}()},
\code{\link{engine_xgboost}()}
}
\concept{engine}

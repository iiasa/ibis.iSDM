% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/engine_gdb.R
\name{engine_gdb}
\alias{engine_gdb}
\title{Use of Gradient Descent Boosting for model estimation}
\usage{
engine_gdb(
  x,
  iter = 2000,
  learning_rate = 0.1,
  empirical_risk = "inbag",
  type = "response",
  ...
)
}
\arguments{
\item{x}{\code{\link[=distribution]{distribution()}} (i.e. \code{\linkS4class{BiodiversityDistribution}}) object.}

\item{iter}{An \code{\link{integer}} giving the number of boosting iterations (Default: \code{2e3L}).}

\item{learning_rate}{A bounded \code{\link{numeric}} value between \code{0} and \code{1}
defining the shrinkage parameter.}

\item{empirical_risk}{method for empirical risk calculation. Available options
are \code{'inbag'}, \code{'oobag'} and \code{'none'}. (Default: \code{'inbag'}).}

\item{type}{The mode used for creating posterior predictions. Either making
\code{"link"}, \code{"response"} or \code{"class"} (Default: \code{"response"}).}

\item{...}{Other variables or control parameters}
}
\value{
An engine.
}
\description{
Gradient descent boosting is an efficient way to optimize any
loss function of a generalized linear or additive model (such as the GAMs
available through the \code{"mgcv"} R-package). It furthermore automatically
regularizes the fit, thus the resulting model only contains the covariates
whose baselearners have some influence on the response. Depending on the type
of the \code{add_biodiversity} data, either poisson process models or
logistic regressions are estimated. If the \code{"only_linear"} term in
\link{train} is set to \code{FALSE}, splines are added to the estimation, thus
providing a non-linear additive inference.
}
\details{
: This package requires the \code{"mboost"} R-package to be
installed. It is in philosophy somewhat related to the \link{engine_xgboost} and
\code{"XGBoost"} R-package, however providing some additional desirable
features that make estimation quicker and particularly useful for spatial
projections. Such as for instance the ability to specifically add spatial
baselearners via \link{add_latent_spatial} or the specification of monotonically
constrained priors via \link{GDBPrior}.
}
\examples{
\dontrun{
# Add GDB as an engine
x <- distribution(background) |> engine_gdb(iter = 1000)
}

}
\references{
\itemize{
\item Hofner, B., Mayr, A., Robinzonov, N., & Schmid, M. (2014). Model-based boosting
in R: a hands-on tutorial using the R package mboost. Computational statistics, 29(1-2), 3-35.
\item Hofner, B., Müller, J., Hothorn, T., (2011). Monotonicity-constrained species
distribution models. Ecology 92, 1895–901.
\item Mayr, A., Hofner, B. and Schmid, M. (2012). The importance of knowing when
to stop - a sequential stopping rule for component-wise gradient boosting.
Methods of Information in Medicine, 51, 178–186.
}
}
\seealso{
Other engine: 
\code{\link{engine_bart}()},
\code{\link{engine_breg}()},
\code{\link{engine_glm}()},
\code{\link{engine_glmnet}()},
\code{\link{engine_inla}()},
\code{\link{engine_inlabru}()},
\code{\link{engine_scampr}()},
\code{\link{engine_stan}()},
\code{\link{engine_xgboost}()}
}
\concept{engine}

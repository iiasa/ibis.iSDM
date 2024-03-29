% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-stan.R
\name{posterior_predict_stanfit}
\alias{posterior_predict_stanfit}
\title{Create a posterior prediction from a rstanfit object}
\usage{
posterior_predict_stanfit(
  obj,
  form,
  newdata,
  mode = "predictor",
  family = NULL,
  offset = NULL,
  draws = NULL
)
}
\arguments{
\item{obj}{A \code{"stanfit"} object (as used by rstan).}

\item{form}{A \code{\link{formula}} object created for the \link{DistributionModel}.}

\item{newdata}{A \link{data.frame} with new data to be used for prediction.}

\item{mode}{A \code{\link{character}} of whether the linear \code{predictor} or the \code{response}
is to be summarized.}

\item{family}{A \code{\link{character}} giving the family for simulating linear response
values (Default: \code{NULL})}

\item{offset}{A \link{vector} with an optionally specified offset.}

\item{draws}{\link{numeric} indicating whether a specific number of draws should be taken.}
}
\description{
This function does simulates from the posterior of a created
stan model, therefore providing a fast and efficient way to project coefficients
obtained from Bayesian models to new/novel contexts.
}
\references{
\itemize{
\item \url{https://medium.com/@alex.pavlakis/making-predictions-from-stan-models-in-r-3e349dfac1ed}.
\item The brms R-package.
}
}

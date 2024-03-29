% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-stan.R
\name{stancode}
\alias{stancode}
\alias{stancode.DistributionModel}
\title{Show the stan code from a trained model}
\usage{
stancode(obj, ...)

stancode.DistributionModel(obj, ...)
}
\arguments{
\item{obj}{Any prepared object.}

\item{...}{not used.}
}
\value{
None.
}
\description{
This helper function shows the code from a trained
\link{DistributionModel} using the \code{\link{engine_stan}}. This function is emulated
after a similar functionality in the brms R-package.
\strong{It only works with models inferred with stan!}
}
\seealso{
rstan, cmdstanr, brms
}
\keyword{engine}

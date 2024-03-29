% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ensemble.R
\name{ensemble_spartial}
\alias{ensemble_spartial}
\alias{ensemble_spartial,ANY-method}
\title{Function to create an ensemble of spartial effects from multiple models}
\usage{
ensemble_spartial(
  ...,
  x.var,
  method = "mean",
  layer = "mean",
  newdata = NULL,
  min.value = NULL,
  normalize = TRUE
)

\S4method{ensemble_spartial}{ANY}(
  ...,
  x.var,
  method = "mean",
  layer = "mean",
  newdata = NULL,
  min.value = NULL,
  normalize = TRUE
)
}
\arguments{
\item{...}{Provided \code{\link{DistributionModel}} objects from which partial responses
can be called. In the future provided data.frames might be supported as well.}

\item{x.var}{A \code{\link{character}} of the variable from which an ensemble is to be
created.}

\item{method}{Approach on how the ensemble is to be created. See details for
options (Default: \code{'mean'}).}

\item{layer}{A \code{\link{character}} of the layer to be taken from each prediction
(Default: \code{'mean'}). If set to \code{NULL} ignore any of the layer names
in ensembles of \code{SpatRaster} objects.}

\item{newdata}{A optional \code{\link{data.frame}} or \code{\link{SpatRaster}} object supplied to
the model (DefaultL \code{NULL}). This object needs to have identical names as the original predictors.}

\item{min.value}{A optional \code{\link{numeric}} stating a minimum value that needs
to be surpassed in each layer before calculating and ensemble
(Default: \code{NULL}).}

\item{normalize}{\code{\link{logical}} on whether the inputs of the ensemble should be
normalized to a scale of 0-1 (Default: \code{TRUE}).}
}
\value{
A \link{SpatRaster} object with the combined partial effects of the supplied models.
}
\description{
Similar to the \code{ensemble()} function, this function creates an
ensemble of partial responses of provided distribution models fitted with
the \code{\link{ibis.iSDM-package}}. Through the \code{layer} parameter it can be
specified which part of the partial prediction should be averaged in an
ensemble (if given). This can be for instance the \emph{mean} prediction and/or
the standard deviation \emph{sd}. Ensemble partial is also being called if more
than one input \code{\link{DistributionModel}} object is provided to \code{partial}.

By default the ensemble of partial responses is created as average across
all models with the uncertainty being the standard deviation of responses.
}
\details{
Possible options for creating an ensemble includes:
\itemize{
\item \code{'mean'} - Calculates the mean of several predictions.
\item \code{'median'} - Calculates the median of several predictions.
}
}
\note{
If a list is supplied, then it is assumed that each entry in the list
is a fitted \code{\link{DistributionModel}} object. Take care not to create an ensemble
of models constructed with different link functions, e.g. logistic vs \link{log}.
By default the response functions of each model are normalized.
}
\examples{
\dontrun{
 # Assumes previously computed models
 ex <- ensemble_spartial(mod1, mod2, mod3, method = "mean")
}

}
\keyword{partial}
\keyword{train}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudoabsence.R
\name{add_pseudoabsence}
\alias{add_pseudoabsence}
\title{Add pseudo-absence points to a point data set}
\usage{
add_pseudoabsence(
  df,
  field_occurrence = "observed",
  template = NULL,
  settings = getOption("ibis.pseudoabsence")
)
}
\arguments{
\item{df}{A \code{\link{sf}}, \code{\link{data.frame}} or \code{\link{tibble}} object containing point data.}

\item{field_occurrence}{A \code{\link{character}} name of the column containing the presence
information (Default: \code{observed}).}

\item{template}{A \code{\link{SpatRaster}} object that is aligned with the predictors
(Default: \code{NULL}). If set to \code{NULL}, then \code{background} in
the \code{\link[=pseudoabs_settings]{pseudoabs_settings()}} has to be a \code{\link{SpatRaster}} object.}

\item{settings}{A \code{\link[=pseudoabs_settings]{pseudoabs_settings()}} objects. Absence settings are
taken from \link{ibis_options} otherwise (Default).}
}
\value{
A \code{\link{data.frame}} containing the newly created pseudo absence points.
}
\description{
For most engines, background or pseudo-absence points are
necessary. The distinction lies in how the absence data are handled. For
\code{\link{poisson}} distributed responses, absence points are considered background
points over which the intensity of sampling (\code{lambda}) is integrated (in
a classical Poisson point-process model).

In contrast in \code{\link{binomial}} distributed responses, the absence information is
assumed to be an adequate representation of the true absences and treated by
the model as such... Here it is advised to specify absence points in a way
that they represent potential true absence, such as for example through
targeted background sampling or by sampling them within/outside a given
range.
}
\details{
A \code{\link[=pseudoabs_settings]{pseudoabs_settings()}} object can be added to setup how absence
points should be sampled. A \code{bias} parameter can be set to specify a
bias layer to sample from, for instance a layer of accessibility. Note that
when modelling several datasets, it might make sense to check across all
datasets whether certain areas are truly absent. By default, the
pseudo-absence points are not sampled in areas in which there are already
presence points.
}
\note{
This method removes all columns from the input \code{df} object other
than the \code{field_occurrence} column and the coordinate columns (which
will be created if not already present).
}
\references{
\itemize{
\item Stolar, J., & Nielsen, S. E. (2015). Accounting for spatially biased sampling effort in
presence‐only species distribution modelling. Diversity and Distributions,
21(5), 595-608.
\item Bird, T.J., Bates, A.E., Lefcheck, J.S., Hill, N.A., Thomson, R.J., Edgar, G.J.,
Stuart-Smith, R.D., Wotherspoon, S., Krkosek, M., Stuart-Smith, J.F. and Pecl, G.T., 2014.
Statistical solutions for error and bias in global citizen science datasets.
Biological Conservation, 173, pp.144-154.
}
}
\keyword{train}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_output.R
\name{load_model}
\alias{load_model}
\alias{load_model,character-method}
\title{Load a pre-computed model}
\usage{
load_model(fname, verbose = getOption("ibis.setupmessages", default = TRUE))

\S4method{load_model}{character}(fname, verbose = getOption("ibis.setupmessages", default = TRUE))
}
\arguments{
\item{fname}{A \code{\link{character}} depicting an output filename.}

\item{verbose}{\code{\link{logical}} indicating whether messages should be shown.
Overwrites \code{getOption("ibis.setupmessages")} (Default: \code{TRUE}).}
}
\value{
A \code{\link{DistributionModel}} object.
}
\description{
The \code{load_model} function (opposed to the \code{write_model}) loads
previous saved \code{\link{DistributionModel}}. It is essentially a wrapper to
\code{\link{readRDS}}.

When models are loaded, they are briefly checked for their validity and
presence of necessary components.
}
\examples{
\dontrun{
# Load model
mod <- load_model("testmodel.rds")

summary(mod)
}

}
\seealso{
write_model
}
\keyword{utils}

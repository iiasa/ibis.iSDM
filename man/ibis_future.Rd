% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{ibis_future}
\alias{ibis_future}
\title{Internal function to enable  (a)synchronous parallel processing}
\usage{
ibis_future(
  plan_exists = FALSE,
  cores = getOption("ibis.nthread", default = 2),
  strategy = getOption("ibis.futurestrategy"),
  workers = NULL
)
}
\arguments{
\item{plan_exists}{A \code{\link{logical}} check on whether an existing \code{\link{future}} plan exists (Default: \code{FALSE}).}

\item{cores}{A \code{\link{numeric}} number stating the number of cores to use.}

\item{strategy}{A \code{\link{character}} denoting the strategy to be used for future.
See help of \code{\link{future}} for options. (Default: \code{"multisession"}).}

\item{workers}{An optional list of remote machines or workers, e.g. \code{"c(remote.server.org)"}.
Alternatively a \code{"cluster"} object can be provided.}
}
\value{
Invisible
}
\description{
This function checks if parallel processing can be set up and enables it.
\strong{Ideally this is done by the user for more control!}
In the package parallelization is usually only used for predictions and projections,
but not for inference in which case parallel inference should be handled by the engine.
}
\details{
Currently supported strategies are:
\itemize{
\item \code{"sequential"} = Resolves futures sequentially in the current R process (Package default).
\item \code{"multisession"} = Resolves futures asynchronously across \code{'cores'} sessions.
\item \code{"multicore"} = Resolves futures asynchronously across on forked processes. Only works on UNIX systems!
\item \code{"cluster"} = Resolves futures asynchronously in sessions on this or more machines.
\item \code{"slurm"} = To be implemented: Slurm linkage via batchtools.
}
}
\note{
The \code{'plan'} set by \link{future} exists after the function has been executed.

If the aim is to parallize across many species, this is better done in a scripted solution.
Make sure not to parallize predictions within existing clusters to avoid out-of-memory
issues.
}
\examples{
\dontrun{
# Starts future job. F in this case is a prediction function.
ibis_future(cores = 4, strategy = "multisession")
}

}
\seealso{
\link{future}
}
\keyword{misc}

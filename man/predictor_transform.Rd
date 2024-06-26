% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-predictors.R
\name{predictor_transform}
\alias{predictor_transform}
\title{Spatial adjustment of environmental predictors and raster stacks}
\usage{
predictor_transform(
  env,
  option,
  windsor_props = c(0.05, 0.95),
  pca.var = 0.8,
  state = NULL,
  method = NULL,
  ...
)
}
\arguments{
\item{env}{A \code{\link{SpatRaster}} or \code{\link{stars}} object.}

\item{option}{A \code{\link{vector}} stating whether predictors should be preprocessed
in any way (Options: \code{'none'}, \code{'scale'}, \code{'norm'}, \code{'windsor'},
\code{'windsor_thresh'}, \code{'percentile'} \code{'pca'}, \code{'revjack'}). See Details.}

\item{windsor_props}{A \code{\link{numeric}} vector specifying the proportions to be
clipped for windsorization (Default: \code{c(.05,.95)}).}

\item{pca.var}{A \code{\link{numeric}} value between \code{>0} and \code{1} stating the
minimum amount of variance to be covered (Default: \code{0.8}).}

\item{state}{A \code{\link{matrix}} with one value per variable (column) providing either a ( \code{stats::mean()}, \code{stats::sd()} )
for each variable in \code{env} for option \code{'scale'} or a range of minimum and maximum values for
option \code{'norm'}. Effectively applies their value range for rescaling. (Default: \code{NULL}).}

\item{method}{As \code{'option'} for more intuitive method setting. Can be left
empty (in this case option has to be set).}

\item{...}{other options (Non specified).}
}
\value{
Returns a adjusted \code{\link{SpatRaster}} object of identical resolution.
}
\description{
This function allows the transformation of provided environmental
predictors (in \code{\link{SpatRaster}} format). A common use case is for instance the
standardization (or scaling) of all predictors prior to model fitting. This
function works both with \code{\link{SpatRaster}} as well as with \code{\link{stars}} objects.
}
\details{
Available options are:
\itemize{
\item \code{'none'} The original layer(s) are returned.
\item \code{'scale'} This run the \code{\link[=scale]{scale()}} function with default settings
(1 Standard deviation) across all predictors. A sensible default to for most model fitting.
\item \code{'norm'} This normalizes all predictors to a range from \code{0-1}.
\item \code{'windsor'} This applies a 'windsorization' to an existing raster layer
by setting the lowest, respectively largest values to the value at a certain
percentage level (e.g. 95\%). Those can be set via the parameter \code{"windsor_props"}.
\item \code{'windsor_thresh'} Same as option 'windsor', however in this case values
are clamped to a thresholds rather than certain percentages calculated on the data.
\item \code{'percentile'} This converts and bins all values into percentiles, e.g.
the top 10\% or lowest 10\% of values and so on.
\item \code{'pca'} This option runs a principal component decomposition of all
predictors (via \code{\link[=prcomp]{prcomp()}}). It returns new predictors resembling all components
in order of the most important ones. Can be useful to reduce collinearity, however
note that this changes all predictor names to 'PCX', where X is the number of
the component. The parameter \code{'pca.var'} can be modified to specify the
minimum variance to be covered by the axes.
\item \code{'revjack'} Removes outliers from the supplied stack via a reverse jackknife
procedure. Identified outliers are by default set to \code{NA}.
}
}
\note{
If future covariates are rescaled or normalized, it is highly recommended to use the
statistical moments on which the models were trained for any variable transformations,
also to ensure that variable ranges are consistent among relative values.
}
\examples{
\dontrun{
# Where x is a SpatRaster
new_x <- predictor_transform(x, option = 'scale')
}

}
\seealso{
predictor_derivate
}
\keyword{utils}

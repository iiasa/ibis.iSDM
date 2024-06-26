% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/threshold.R
\name{threshold}
\alias{threshold}
\alias{threshold,ANY-method}
\alias{threshold,SpatRaster-method}
\alias{threshold,BiodiversityScenario-method}
\title{Threshold a continuous prediction to a categorical layer}
\usage{
threshold(
  obj,
  method = "mtp",
  value = NULL,
  point = NULL,
  field_occurrence = "observed",
  format = "binary",
  return_threshold = FALSE,
  ...
)

\S4method{threshold}{ANY}(
  obj,
  method = "mtp",
  value = NULL,
  point = NULL,
  field_occurrence = "observed",
  format = "binary",
  return_threshold = FALSE,
  ...
)

\S4method{threshold}{SpatRaster}(
  obj,
  method = "fixed",
  value = NULL,
  point = NULL,
  field_occurrence = "observed",
  format = "binary",
  return_threshold = FALSE
)

\S4method{threshold}{BiodiversityScenario}(
  obj,
  method = "mtp",
  value = NULL,
  point = NULL,
  field_occurrence = "observed",
  format = "binary",
  return_threshold = FALSE,
  ...
)
}
\arguments{
\item{obj}{A \link{BiodiversityScenario} object to which an existing threshold is
to be added.}

\item{method}{A specifc method for thresholding. See details for available options.}

\item{value}{A \code{\link{numeric}} value specifying the specific threshold for scenarios
(Default: \code{NULL} Grab from object).}

\item{point}{A \code{\link{sf}} object containing observational data used for model training.}

\item{field_occurrence}{A \code{\link{character}} location of biodiversity point records.}

\item{format}{\code{\link{character}} indication of whether \code{"binary"}, \code{"normalize"}
or \code{"percentile"} formatted thresholds are to be created (Default: \code{"binary"}).
Also see Muscatello et al. (2021).}

\item{return_threshold}{Should threshold value be returned instead (Default: \code{FALSE})}

\item{...}{Any other parameter. Used to fetch value if set somehow.}
}
\value{
A \link{SpatRaster} if a \link{SpatRaster} object as input. Otherwise the threshold
is added to the respective \code{\link{DistributionModel}} or \code{\link{BiodiversityScenario}} object.
}
\description{
It is common in many applications of species distribution
modelling that estimated continuous suitability surfaces are converted into
discrete representations of where suitable habitat might or might not exist.
This so called \emph{threshold'ing} can be done in various ways which are further
described in the details.

In case a \code{\link{SpatRaster}} is provided as input in this function for
\code{obj}, it is furthermore necessary to provide a \code{\link{sf}} object for
validation as there is no \code{\link{DistributionModel}} to read this information
from.

\strong{Note:} This of course also allows to estimate the threshold based on withheld data, for instance
those created from an a-priori cross-validation procedure.

For \code{\link{BiodiversityScenario}} objects, adding this function to the processing
pipeline stores a threshold attribute in the created \link{scenario} object.

For \link{BiodiversityScenario} objects a set \code{threshold()} simply indicates that
the projection should create and use thresholds as part of the results.
The threshold values for this are either taken from the provided model or
through an optional provide parameter \code{value}.

If instead the aim is to apply thresholds to each step of the suitability
projection, see \code{\link[=add_constraint_threshold]{add_constraint_threshold()}}.
}
\details{
The following options are currently implemented:
\itemize{
\item \code{'fixed'} = applies a single pre-determined threshold. Requires \code{value}
to be set.
\item \code{'mtp'} = minimum training presence is used to find and set the lowest
predicted suitability for any occurrence point.
\item \code{'percentile'} = For a percentile threshold. A \code{value} as parameter
has to be set here.
\item \code{'min.cv'} = Threshold the raster so to minimize the coefficient of
variation (cv) of the posterior. Uses the lowest tercile of the cv in space.
Only feasible with Bayesian engines.
\item \code{'TSS'} = Determines the optimal TSS (True Skill Statistic). Requires
the \code{"modEvA"} package to be installed.
\item \code{'kappa'} = Determines the optimal kappa value (Kappa). Requires the
\code{"modEvA"} package to be installed.
\item \code{'F1score'} = Determines the optimal F1score (also known as Sorensen
similarity). Requires the \code{"modEvA"} package to be installed.
\item \code{'F1score'} = Determines the optimal sensitivity of presence records.
Requires the \code{"modEvA"} package to be installed.
\item \code{'Sensitivity'} = Determines the optimal sensitivity of presence records.
Requires the \code{"modEvA"} package to be installed.
\item \code{'Specificity'} = Determines the optimal sensitivity of presence records.
Requires the \code{"modEvA"} package to be installed.
\item \code{'AUC'} = Determines the optimal AUC of presence records. Requires the
\code{"modEvA"} package to be installed.
}
}
\examples{
\dontrun{
 # Where mod is an estimated DistributionModel
 tr <- threshold(mod)
 tr$plot_threshold()
}

}
\references{
\itemize{
\item Lawson, C.R., Hodgson, J.A., Wilson, R.J., Richards, S.A., 2014. Prevalence,
thresholds and the performance of presence-absence models. Methods Ecol. Evol.
5, 54–64. https://doi.org/10.1111/2041-210X.12123
\item Liu, C., White, M., Newell, G., 2013. Selecting thresholds for the prediction
of species occurrence with presence-only data. J. Biogeogr. 40, 778–789. https://doi.org/10.1111/jbi.12058
\item Muscatello, A., Elith, J., Kujala, H., 2021. How decisions about fitting
species distribution models affect conservation outcomes. Conserv. Biol. 35, 1309–1320.
https://doi.org/10.1111/cobi.13669
}
}
\seealso{
\code{"modEvA"}
}

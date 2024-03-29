% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{distribution}
\alias{distribution}
\alias{distribution,SpatRaster-method}
\alias{distribution,sf-method}
\title{Create distribution modelling procedure}
\usage{
distribution(
  background,
  limits = NULL,
  limits_method = "none",
  mcp_buffer = 0,
  limits_clip = FALSE
)

\S4method{distribution}{SpatRaster}(
  background,
  limits = NULL,
  limits_method = "none",
  mcp_buffer = 0,
  limits_clip = FALSE
)

\S4method{distribution}{sf}(
  background,
  limits = NULL,
  limits_method = "none",
  mcp_buffer = 0,
  limits_clip = FALSE
)
}
\arguments{
\item{background}{Specification of the modelling background. Must be a
\code{\link{SpatRaster}} or \code{\link{sf}} object.}

\item{limits}{A \code{\link{SpatRaster}} or \code{\link{sf}} object that limits the prediction
surface when intersected with input data (Default: \code{NULL}).}

\item{limits_method}{A \code{\link{character}} of the method used for hard limiting a
projection. Available options are \code{"none"} (Default), \code{"zones"}
or \code{"mcp"}.}

\item{mcp_buffer}{A \code{\link{numeric}} distance to buffer the mcp (Default \code{0}).
Only used if \code{"mcp"} is used.}

\item{limits_clip}{\code{\link{logical}} Should the limits clip all predictors before
fitting a model (\code{TRUE}) or just the prediction (\code{FALSE}, default).}
}
\value{
\code{\linkS4class{BiodiversityDistribution}} object containing data for building
a biodiversity distribution modelling problem.
}
\description{
This function creates an object that contains all the data,
parameters and settings for building an (integrated) species distribution
model. Key functions to add data are \code{\link{add_biodiversity_poipo}} and the like,
\code{\link{add_predictors}}, \code{\link{add_latent_spatial}}, \code{\link{engine_glmnet}} or similar,
\code{\link{add_priors}} and \code{\link{add_offset}}. It creates a prototype
\code{\link{BiodiversityDistribution}} object with its own functions. After setting
input data and parameters, model predictions can then be created via the
\link{train} function and predictions be created.

Additionally, it is possible to specify a \code{"limit"} to any predictions
conducted on the background. This can be for instance a buffered layer by a
certain dispersal distance (Cooper and Soberon, 2018) or a categorical layer
representing biomes or soil conditions. Another option is to create a
constraint by constructing a minimum convex polygon (MCP) using the supplied
biodiversity data. This option can be enabled by setting
\code{"limits_method"} to \code{"mcp"}. It is also possible to provide a
small buffer to constructed MCP that way. See the frequently asked question
(FAQ) section on the homepage for more information.

See \strong{Details} for a description of the internal functions available to
modify or summarize data within the created object.

\strong{Note that any model requires at minimum a single added biodiversity dataset
as well as a specified engine.}
}
\details{
This function creates a \code{\linkS4class{BiodiversityDistribution}} object
that in itself contains other functions and stores parameters and
(pre-)processed data. A full list of functions available can be queried via
\code{"names(object)"}. Some of the functions are not intended to be
manipulated directly, but rather through convenience functions (e.g.
\code{"object$set_predictors()"}). Similarly other objects are stored in the
\code{\linkS4class{BiodiversityDistribution}} object that have their own functions as
well and can be queried (e.g. \code{"names(object)"}). For a list of
functions see the reference documentation. By default, if some datasets are
not set, then a \code{"Waiver"} object is returned instead.

The following objects can be stored:
\itemize{
\item \code{object$biodiversity} A \code{\link{BiodiversityDatasetCollection}} object with
the added biodiversity data.
\item \code{object$engine} An \code{"engine"} object (e.g. \code{\link[=engine_inlabru]{engine_inlabru()}})
with function depended on the added engine.
\item \code{object$predictors} A \code{\link{PredictorDataset}} object with all set predictions.
\item \code{object$priors} A \code{\link{PriorList}} object with all specified priors.
\item \code{object$log} A \code{\link{Log}} object that captures.
}

Useful high-level functions to address those objects are for instance:
\itemize{
\item \code{object$show()} A generic summary of the \code{\linkS4class{BiodiversityDistribution}}
object contents. Can also be called via \link{print}.
\item \code{object$get_biodiversity_equations()} Lists the equations used for each
biodiversity dataset with given id. Defaults to all predictors.
\item \code{object$get_biodiversity_types()} Lists the type of each specified
biodiversity dataset with given id.
\item \code{object$get_extent()} Outputs the \link[terra:ext]{terra::ext} of the modelling region.
\item \code{object$show_background_info()} Returns a \code{\link{list}} with the \link[terra:ext]{terra::ext}
and the \link[terra:crs]{terra::crs}.
\item \code{object$get_extent_dimensions()} Outputs the \link[terra:ext]{terra::ext} dimension by
calling the \code{"extent_dimensions()"} function.
\item \code{object$get_predictor_names()} Returns a \link{character} vector with the
names of all added predictors.
\item \code{object$get_prior_variables()} Returns a description of \code{\link{priors}} added.
}

There are other functions as well but those are better accessed through their
respective wrapper functions.
}
\examples{
# Load background raster
background <- terra::rast(system.file("extdata/europegrid_50km.tif",package = "ibis.iSDM"))
# Define model
x <- distribution(background)
x

}
\references{
\itemize{
\item Fletcher, R.J., Hefley, T.J., Robertson, E.P., Zuckerberg, B., McCleery, R.A.,
Dorazio, R.M., (2019) A practical guide for combining data to model species
distributions. Ecology 100, e02710. https://doi.org/10.1002/ecy.2710
\item Cooper, Jacob C., and Jorge Soberón. "Creating individual accessible area
hypotheses improves stacked species distribution model performance." Global
Ecology and Biogeography 27, no. 1 (2018): 156-165.
}
}
\seealso{
\code{\link{BiodiversityDistribution}} and other classes.
}

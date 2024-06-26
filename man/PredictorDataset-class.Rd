% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class-predictors.R
\name{PredictorDataset-class}
\alias{PredictorDataset-class}
\alias{PredictorDataset}
\title{PredictorDataset class description}
\description{
This class describes the PredictorDataset and is used to store covariates within.
}
\seealso{
\code{\link[=predictor_derivate]{predictor_derivate()}}

\code{\link[=predictor_transform]{predictor_transform()}}

\code{\link[=predictor_transform]{predictor_transform()}}
}
\keyword{classes}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{id}}{The id for this collection as \code{\link{character}}.}

\item{\code{data}}{A predictor dataset usually as \code{\link{SpatRaster}}.}

\item{\code{name}}{A name for this object.}

\item{\code{transformed}}{Saves whether the predictors have been transformed somehow.}

\item{\code{timeperiod}}{A timeperiod field}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-PredictorDataset-new}{\code{PredictorDataset$new()}}
\item \href{#method-PredictorDataset-print}{\code{PredictorDataset$print()}}
\item \href{#method-PredictorDataset-get_name}{\code{PredictorDataset$get_name()}}
\item \href{#method-PredictorDataset-get_id}{\code{PredictorDataset$get_id()}}
\item \href{#method-PredictorDataset-get_names}{\code{PredictorDataset$get_names()}}
\item \href{#method-PredictorDataset-get_predictor_names}{\code{PredictorDataset$get_predictor_names()}}
\item \href{#method-PredictorDataset-get_data}{\code{PredictorDataset$get_data()}}
\item \href{#method-PredictorDataset-get_time}{\code{PredictorDataset$get_time()}}
\item \href{#method-PredictorDataset-get_projection}{\code{PredictorDataset$get_projection()}}
\item \href{#method-PredictorDataset-get_resolution}{\code{PredictorDataset$get_resolution()}}
\item \href{#method-PredictorDataset-crop_data}{\code{PredictorDataset$crop_data()}}
\item \href{#method-PredictorDataset-mask}{\code{PredictorDataset$mask()}}
\item \href{#method-PredictorDataset-set_data}{\code{PredictorDataset$set_data()}}
\item \href{#method-PredictorDataset-rm_data}{\code{PredictorDataset$rm_data()}}
\item \href{#method-PredictorDataset-show}{\code{PredictorDataset$show()}}
\item \href{#method-PredictorDataset-summary}{\code{PredictorDataset$summary()}}
\item \href{#method-PredictorDataset-has_derivates}{\code{PredictorDataset$has_derivates()}}
\item \href{#method-PredictorDataset-is_transformed}{\code{PredictorDataset$is_transformed()}}
\item \href{#method-PredictorDataset-get_transformed_params}{\code{PredictorDataset$get_transformed_params()}}
\item \href{#method-PredictorDataset-length}{\code{PredictorDataset$length()}}
\item \href{#method-PredictorDataset-ncell}{\code{PredictorDataset$ncell()}}
\item \href{#method-PredictorDataset-plot}{\code{PredictorDataset$plot()}}
\item \href{#method-PredictorDataset-clone}{\code{PredictorDataset$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-new"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-new}{}}}
\subsection{Method \code{new()}}{
Initializes the object and creates an empty list
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$new(id, data, transformed = FALSE, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{id}}{The id for this collection as \code{\link{character}}.}

\item{\code{data}}{A predictor dataset usually as \code{\link{SpatRaster}}.}

\item{\code{transformed}}{A \code{\link{logical}} flag if predictors have been transformed. Assume not.}

\item{\code{...}}{Any other parameters found.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
NULL
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-print"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-print}{}}}
\subsection{Method \code{print()}}{
Print the names and properties of all Biodiversity datasets contained within
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$print(format = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{format}}{A \code{\link{logical}} flag on whether a message should be printed.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A message on screen
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-get_name"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-get_name}{}}}
\subsection{Method \code{get_name()}}{
Return name of this object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$get_name()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Default \code{\link{character}} name.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-get_id"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-get_id}{}}}
\subsection{Method \code{get_id()}}{
Get Id of this object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$get_id()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Default \code{\link{character}} name.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-get_names"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-get_names}{}}}
\subsection{Method \code{get_names()}}{
Get names of data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$get_names()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
\code{\link{character}} names of the data value.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-get_predictor_names"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-get_predictor_names}{}}}
\subsection{Method \code{get_predictor_names()}}{
Alias for get_names
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$get_predictor_names()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
\code{\link{character}} names of the data value.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-get_data"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-get_data}{}}}
\subsection{Method \code{get_data()}}{
Get a specific dataset
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$get_data(df = FALSE, na.rm = TRUE, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{df}}{\code{\link{logical}} on whether data is to be returned as \code{\link{data.frame}}.}

\item{\code{na.rm}}{\code{\link{logical}} if \code{NA} is to be removed from data.frame.}

\item{\code{...}}{Any other parameters passed on.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A \code{\link{SpatRaster}} or \code{\link{data.frame}}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-get_time"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-get_time}{}}}
\subsection{Method \code{get_time()}}{
Get time dimension of object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$get_time(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Any other parameters passed on.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A \code{\link{vector}} with the time dimension of the dataset.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-get_projection"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-get_projection}{}}}
\subsection{Method \code{get_projection()}}{
Get Projection
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$get_projection()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A \code{\link{vector}} with the geographical projection of the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-get_resolution"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-get_resolution}{}}}
\subsection{Method \code{get_resolution()}}{
Get Resolution
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$get_resolution()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A \code{\link{numeric}} \code{\link{vector}} with the spatial resolution of the data.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-crop_data"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-crop_data}{}}}
\subsection{Method \code{crop_data()}}{
Utility function to clip the predictor dataset by another dataset
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$crop_data(pol)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{pol}}{A \code{\link{sf}} object used for cropping the data}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Invisibile TRUE
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-mask"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-mask}{}}}
\subsection{Method \code{mask()}}{
Utility function to mask the predictor dataset by another dataset
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$mask(mask, inverse = FALSE, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{mask}}{A \code{SpatRaster} or \code{sf} object.}

\item{\code{inverse}}{A \code{logical} flag if the inverse should be masked instead.}

\item{\code{...}}{Any other parameters passed on to masking.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Invisible
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-set_data"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-set_data}{}}}
\subsection{Method \code{set_data()}}{
Add a new Predictor dataset to this collection
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$set_data(x, value)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\code{\link{character}} of the new name to be stored.}

\item{\code{value}}{A new \code{\link{SpatRaster}} object.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
This object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-rm_data"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-rm_data}{}}}
\subsection{Method \code{rm_data()}}{
Remove a specific Predictor by name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$rm_data(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\code{\link{character}} of the predictor name to be removed.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Invisible
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-show"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-show}{}}}
\subsection{Method \code{show()}}{
Alias for print method
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$show()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Invisible
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-summary"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-summary}{}}}
\subsection{Method \code{summary()}}{
Collect info statistics with optional decimals
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$summary(digits = 2)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{digits}}{\code{\link{numeric}} Giving the rounding precision}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A \code{\link{data.frame}} summarizing the data.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-has_derivates"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-has_derivates}{}}}
\subsection{Method \code{has_derivates()}}{
Indication if there are any predictors that are derivates of outers
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$has_derivates()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A \code{\link{logical}} flag.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-is_transformed"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-is_transformed}{}}}
\subsection{Method \code{is_transformed()}}{
Predictors have been transformed?
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$is_transformed()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A \code{\link{logical}} flag.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-get_transformed_params"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-get_transformed_params}{}}}
\subsection{Method \code{get_transformed_params()}}{
Get transformation params.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$get_transformed_params()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A \code{\link{matrix}} flag.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-length"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-length}{}}}
\subsection{Method \code{length()}}{
Number of Predictors in object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$length()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A \code{\link{numeric}} estimate
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-ncell"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-ncell}{}}}
\subsection{Method \code{ncell()}}{
Number of cells or values in object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$ncell()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A \code{\link{numeric}} estimate
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-plot"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-plot}{}}}
\subsection{Method \code{plot()}}{
Basic Plotting function
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$plot()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A graphical interpretation of the predictors in this object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PredictorDataset-clone"></a>}}
\if{latex}{\out{\hypertarget{method-PredictorDataset-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PredictorDataset$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}

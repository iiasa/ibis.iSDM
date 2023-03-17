#' @include utils.R
NULL

#' Threshold a continuous prediction to a categorical layer
#'
#' @description
#' It is common in many applications of species distribution modelling that estimated
#' continuous suitability surfaces are converted into discrete representations of where
#' suitable habitat might or might not exist. This so called *threshold'ing*
#' can be done in various ways which are further described in the details.
#'
#' In case a [RasterLayer] or [RasterBrick] is provided as input in this function
#' for \code{obj}, it is furthermore necessary to provide a [`sf`] object for validation as
#' there is no [`DistributionModel`] to read this information from.
#' **Note:** This of course also allows to estimate the threshold based on withheld data, for instance
#' those created from an a-priori cross-validation procedure.
#'
#' For [`BiodiversityScenario`] objects, adding this function to the processing pipeline
#' stores a threshold attribute in the created [scenario] object.
#'
#' @param obj A trained [`DistributionModel`] or alternatively a [`Raster`] object.
#' @param method A specifc method for thresholding. See details for available options.
#' @param value A [`numeric`] value for thresholding if method is fixed (Default: \code{NULL}).
#' @param poi A [`sf`] object containing observational data used for model training.
#' @param format [`character`] indication of whether \code{"binary"}, \code{"normalize"} or \code{"percentile"}
#' formatted thresholds are to be created (Default: \code{"binary"}). Also see Muscatello et al. (2021).
#' @param ... other parameters not yet set.
#' @param return_threshold Should threshold value be returned instead (Default: \code{FALSE})
#' @details
#' The following options are currently implemented:
#' * \code{'fixed'} = applies a single pre-determined threshold. Requires \code{value} to be set.
#' * \code{'mtp'} = minimum training presence is used to find and set the lowest predicted suitability for any occurrence point.
#' * \code{'percentile'} = For a percentile threshold. A \code{value} as parameter has to be set here.
#' * \code{'min.cv'} = Threshold the raster so to minimize the coefficient of variation (cv) of the posterior. Uses the lowest tercile of the cv in space. Only feasible with Bayesian engines.
#' * \code{'TSS'} = Determines the optimal TSS (True Skill Statistic). Requires the [modEvA] package to be installed.
#' * \code{'kappa'} = Determines the optimal kappa value (Kappa). Requires the [modEvA] package to be installed.
#' * \code{'F1score'} = Determines the optimal F1score (also known as Sorensen similarity). Requires the [modEvA] package to be installed.
#' * \code{'F1score'} = Determines the optimal sensitivity of presence records. Requires the [modEvA] package to be installed.
#' * \code{'Sensitivity'} = Determines the optimal sensitivity of presence records. Requires the [modEvA] package to be installed.
#' * \code{'Specificity'} = Determines the optimal sensitivity of presence records. Requires the [modEvA] package to be installed.
#' @name threshold
#' @references
#' * Lawson, C.R., Hodgson, J.A., Wilson, R.J., Richards, S.A., 2014. Prevalence, thresholds and the performance of presence-absence models. Methods Ecol. Evol. 5, 54–64. https://doi.org/10.1111/2041-210X.12123
#' * Liu, C., White, M., Newell, G., 2013. Selecting thresholds for the prediction of species occurrence with presence-only data. J. Biogeogr. 40, 778–789. https://doi.org/10.1111/jbi.12058
#' * Muscatello, A., Elith, J., Kujala, H., 2021. How decisions about fitting species distribution models affect conservation outcomes. Conserv. Biol. 35, 1309–1320. https://doi.org/10.1111/cobi.13669
#' @seealso [modEvA]
#' @returns A [RasterLayer] if used with a [Raster] object as input.
#' Otherwise the threshold is added to the respective [`DistributionModel`] or [`BiodiversityScenario`] object.
#' @examples
#' \dontrun{
#'  # Where mod is an estimated DistributionModel
#'  tr <- threshold(mod)
#'  tr$plot_threshold()
#' }
#' @export
NULL

#' @name threshold
#' @rdname threshold
#' @exportMethod threshold
#' @export
methods::setGeneric(
  "threshold",
  signature = methods::signature("obj", "method", "value"),
  function(obj, method = 'mtp', value = NULL, poi = NULL,  format = "binary", return_threshold = FALSE, ...) standardGeneric("threshold"))

#' Generic threshold with supplied DistributionModel object
#' @name threshold
#' @rdname threshold
#' @usage \S4method{threshold}{ANY}(obj)
methods::setMethod(
  "threshold",
  methods::signature(obj = "ANY"),
  function(obj, method = 'mtp', value = NULL, format = "binary", return_threshold = FALSE, ...) {
    assertthat::assert_that(any( class(obj) %in% getOption('ibis.engines') ),
                            is.character(method),
                            is.null(value) || is.numeric(value),
                            is.character(format)
    )
    # Check other and add legacy handling
    dots <- list(...)
    if("truncate" %in% names(dots)) format <- ifelse(dots[[truncate]],"normalize", "binary")
    format <- match.arg(format, c("binary", "normalize", "percentile"), several.ok = FALSE)

    # Get prediction raster
    ras <- obj$get_data('prediction')
    # Get model object
    model <- obj$model

    # Check that the object actually contains a prediction
    assertthat::assert_that(
      is.Raster(ras),
      !is.Waiver(ras),
      msg = 'No fitted prediction in object!'
    )
    # Matching for correct method
    method <- match.arg(method, c('fixed','mtp','percentile','min.cv',
                                           'TSS','kappa','F1score','Sensitivity','Specificity'), several.ok = FALSE)

    # If method is min.cv, check that posterior is accessible
    if(method == "min.cv") assertthat::assert_that("cv" %in% names(ras), msg = "Method min.cv requires a posterior prediction and coefficient of variation!")

    # Get all point data in distribution model
    poi <- do.call(sf:::rbind.sf,
                     lapply(obj$model$biodiversity, function(y){
                       o <- guess_sf(y$observations)
                       o$name <- y$name; o$type <- y$type
                       subset(o, select = c('observed', "name", "type", "geometry"))
                     } )
    ) |> tibble::remove_rownames()
    suppressWarnings(
      poi <- sf::st_set_crs(poi, value = sf::st_crs(obj$get_data('prediction')))
    )

    # If TSS or kappa is chosen, check whether there is poipa data among the sources
    if((!any(poi$observed==0) & method %in% c('TSS','kappa','F1score','Sensitivity','Specificity')) || length(unique(poi$name)) > 1){
      if(getOption('ibis.setupmessages')) myLog('[Threshold]','red','Threshold method needs absence-data. Generating some now...')
      bg <- raster::rasterize(obj$model$background, emptyraster(obj$get_data('prediction')))
      abs <- add_pseudoabsence(df = poi,
                                   field_occurrence = 'observed',
                                   template = bg,
                                   # Assuming that settings are comparable among objects
                                   settings = model$biodiversity[[1]]$pseudoabsence_settings
                                   )

      abs <- subset(abs, select = c('x','y'));abs$observed <- 0
      abs <- guess_sf(abs)
      abs$name <- 'Background point'; abs$type <- "generated"
      suppressWarnings(
        abs <- sf::st_set_crs(abs, value = sf::st_crs(obj$get_data('prediction')))
      )
      poi <- subset(poi, select = c("observed", "name", "type","geometry"))
      abs <- subset(abs, select = c("observed", "name", "type","geometry"))
      poi <- rbind(poi, abs);rm(abs)
    }

    # Convert to sf
    if(!inherits(poi,"sf")){ poi <- guess_sf(poi) }

    # Now self call threshold
    out <- threshold(ras, method = method, value = value, poi = poi, format = format,...)
    assertthat::assert_that(is.Raster(out))
    # Add result to new obj
    new_obj <- obj
    if(inherits(out,'RasterLayer')){
      new_obj <- new_obj$set_data(names(out), out)
    } else if(inherits(out,'RasterStack')) {
      # When stack loop through and add
      new_obj <- new_obj$set_data(paste0("threshold_", method), out)
    }
    # Return altered object
    return(new_obj)
  }
)

#' @noRd
#' @keywords internal
.stackthreshold <- function(obj, method = 'fixed', value = NULL,
                            poi = NULL, format = "binary", return_threshold = FALSE, ...) {
  assertthat::assert_that(is.Raster(obj),
                          is.character(method),
                          inherits(poi,'sf'),
                          is.null(value) || is.numeric(value),
                          is.character(format)
  )
  # Match format
  format <- match.arg(format, c("binary", "normalize", "percentile"), several.ok = FALSE)

  # Apply threshold on each entry
  if(return_threshold){
    # Return the threshold directly
    out <- vector()
    for(i in names(obj)) out <- c(out, threshold(obj[[i]], method = method,
                                                                value = value, poi = poi,  format = format, return_threshold = return_threshold, ...) )
    names(out) <- names(obj)
  } else {
    # Return the raster instead
    out <- raster::stack()
    if(method == "min.cv"){
      # If the coefficient of variation is to be minmized, mask first all values with the threshold only
      assertthat::assert_that(raster::nlayers(obj)>2, "sd" %in% names(obj))
      # Get global coefficient of variation
      errortr <- quantile(obj[["cv"]], .3)
      assertthat::assert_that(is.numeric(errortr))
      # Create mask
      mm <- obj[["cv"]]
      mm[mm > errortr] <- NA
      obj <- raster::mask(obj, mm); rm(mm)
      # Set the value to errortr
      value <- errortr
    }
    # Now loop
    for(i in names(obj)) out <- raster::addLayer(out, threshold(obj[[i]], method = method,
                                                                  value = value, poi = poi, format = format, return_threshold = return_threshold, ...) )
  }
  return(out)
}

#' @name threshold
#' @rdname threshold
#' @inheritParams threshold
#' @usage \S4method{threshold}{RasterBrick}(obj)
methods::setMethod("threshold",methods::signature(obj = "RasterBrick"),.stackthreshold)
#' @usage \S4method{threshold}{RasterStack}(obj)
methods::setMethod("threshold",methods::signature(obj = "RasterStack"),.stackthreshold)

#' @name threshold
#' @rdname threshold
#' @usage \S4method{threshold}{RasterLayer}(obj)
methods::setMethod(
  "threshold",
  methods::signature(obj = "RasterLayer"),
  function(obj, method = 'fixed', value = NULL, poi = NULL, format = "binary", return_threshold = FALSE, plot = FALSE) {
    assertthat::assert_that(is.Raster(obj),
                            inherits(obj,'RasterLayer'),
                            is.character(method),
                            is.null(value) || is.numeric(value),
                            is.character(format)
    )
    # Match format
    format <- match.arg(format, c("binary", "normalize", "percentile"), several.ok = FALSE)

    # If poi is set, try to convert sf
    if(!is.null(poi)) try({poi <- sf::st_as_sf(poi)}, silent = TRUE)
    assertthat::assert_that(is.null(poi) || inherits(poi,'sf'))

    # If observed is a factor, convert to numeric
    if(is.factor(poi$observed)){
      poi$observed <- as.numeric(as.character( poi$observed ))
    }

    # Match to correct spelling mistakes
    method <- match.arg(method, c('fixed','mtp','percentile','min.cv',
                                           'TSS','kappa','F1score','Sensitivity','Specificity'), several.ok = FALSE)

    # Check that raster has at least a mean prediction in name
    if(!is.null(poi)) {
      assertthat::assert_that(unique(sf::st_geometry_type(poi)) %in% c('POINT','MULTIPOINT'))
      assertthat::assert_that(utils::hasName(poi, 'observed'))
      poi_pres <- subset(poi, observed > 0) # Remove any eventual absence data for a poi_pres evaluation
    }
    # Get the raster layer
    raster_thresh <- obj

    # Specify by type:
    if(method == "fixed"){
      # Fixed threshold. Confirm to be set
      assertthat::assert_that(is.numeric(value), msg = 'Fixed value is missing!')
      tr <- value
    } else if(method == "mtp"){
      # minimum training presence
      pointVals <- raster::extract(raster_thresh, poi_pres) # Extract point only estimates
      # Minimum threshold
      tr <- min( stats::na.omit(pointVals) )

    } else if(method == "percentile"){
      pointVals <- raster::extract(raster_thresh, poi_pres) # Extract point only estimates
      pointVals <- subset(pointVals, stats::complete.cases(pointVals)) # Remove any NA or NAN data here
      # percentile training threshold
      if(is.null(value)) value <- 0.1 # If value is not set, use 10%
      if(length(pointVals) < 10) {
        perc <- floor(length(pointVals) * (1 - value))
      } else {
        perc <- ceiling(length(pointVals) * (1 - value))
      }
      tr <- rev(sort(pointVals))[perc] # Percentile threshold

    } else if(method == "min.cv"){
      assertthat::assert_that(!is.null(value),msg = "Global minimum cv needs to be set!")
      pointVals <- raster::extract(raster_thresh, poi_pres) # Extract point only estimates

      # Get standard deviation and calculate percentile
      tr <- min( stats::na.omit(pointVals) )
      names(tr) <- "tr"
      names(value) <- "min.cv"
      # Combine as a vector
      tr <- c(tr, value)

    } else {
      # Optimized threshold statistics using the modEvA package
      # FIXME: Could think of porting these functions but too much effort for now. Rather have users install the package here
      check_package("modEvA")
      # Assure that point data is correctly specified
      assertthat::assert_that(inherits(poi, 'sf'), utils::hasName(poi, 'observed'))
      poi$observed <- ifelse(poi$observed>1,1,poi$observed) # Ensure that observed is <=1
      assertthat::assert_that(all( unique(poi$observed) %in% c(0,1) ))

      # Re-extract point vals but with the full dataset
      pointVals <- raster::extract(raster_thresh, poi)
      assertthat::assert_that(length(pointVals)>2)
      # Calculate the optimal thresholds
      suppressWarnings(
        opt <- modEvA::optiThresh(obs = poi$observed, pred = pointVals,
                                  measures = c("TSS","kappa","F1score","Misclass","Omission","Commission",
                                               "Sensitivity","Specificity"),
                                  optimize = "each", plot = plot)
      )
      if(method %in% opt$optimals.each$measure){
        tr <- opt$optimals.each$threshold[which(opt$optimals.each$measure==method)]
      } else {
        # Returning a collection of them as vector
        tr <- opt$optimals.each$threshold; names(tr) <- opt$optimals.each$measure
      }
    }
    # Security check
    assertthat::assert_that(is.numeric(tr) || is.vector(tr))

    # -- Threshold -- #
    if(return_threshold){
      names(tr) <- method
      return(tr)
    } else {
      # Finally threshold the raster
      # Process depending on format
      if(format == "binary"){
        # Default is to create a binary presence-absence. Otherwise truncated hinge
        raster_thresh[raster_thresh < tr[1]] <- 0
        raster_thresh[raster_thresh >= tr[1]] <- 1
        raster_thresh <- raster::asFactor(raster_thresh)
      } else if(format == "normalize"){
        raster_thresh[raster_thresh < tr[1]] <- NA
        # If truncate, ensure that resulting values are normalized
        raster_thresh <- predictor_transform(raster_thresh, option = "norm")
        raster_thresh[is.na(raster_thresh)] <- 0
        raster_thresh <- raster::mask(raster_thresh, obj)
        base::attr(raster_thresh, 'truncate') <- TRUE

        base::attr(raster_thresh, 'truncate') <- TRUE # Legacy truncate attribute
      } else if(format == "percentile") {
        raster_thresh[raster_thresh < tr[1]] <- NA
        raster_thresh <- predictor_transform(raster_thresh, option = "percentile")
        raster_thresh <- raster::mask(raster_thresh, obj)
        base::attr(raster_thresh, 'truncate') <- TRUE
      }
      names(raster_thresh) <- paste0('threshold_',names(obj),'_',method)
      # Assign attributes
      base::attr(raster_thresh, 'method') <- method
      base::attr(raster_thresh, 'format') <- format
      base::attr(raster_thresh, 'threshold') <- tr
    }
    # Return result
    return(raster_thresh)
  }
)

#### For scenarios ####

#' Thresholds in scenario estimation
#'
#' @name threshold
#' @inheritParams threshold
#' @rdname threshold
#' @usage \S4method{threshold}{BiodiversityScenario}(obj)
methods::setMethod(
  "threshold",
  methods::signature(obj = "BiodiversityScenario"),
  function(obj, tr = new_waiver(), ...) {
    # Assert that predicted raster is present
    assertthat::assert_that( is.Raster(obj$get_model()$get_data('prediction')) )
    # Unless set, check
    if(is.Waiver(tr)){
      # Check that a threshold layer is available and get the methods and data from it
      assertthat::assert_that( length( grep('threshold', obj$get_model()$show_rasters()) ) >0 ,
                               msg = 'Call \' threshold \' for prediction first!')
      # Get threshold layer
      tr_lyr <- grep('threshold', obj$get_model()$show_rasters(),value = TRUE)
      if(length(tr_lyr)>1) warning("There appear to be multiple thresholds. Using the first one.")
      ras_tr <- obj$get_model()$get_data( tr_lyr[1] )
      tr <- attr(ras_tr[[1]], 'threshold')
      names(tr) <- attr(ras_tr[[1]], 'method')
    } else {
      assertthat::assert_that(is.numeric(tr))
    }
    bdproto(NULL, obj, threshold = tr)
  }
)

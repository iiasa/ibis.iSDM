#' Threshold a continuous prediction to a categorical layer
#'
#' @description It is common in many applications of species distribution
#' modelling that estimated continuous suitability surfaces are converted into
#' discrete representations of where suitable habitat might or might not exist.
#' This so called *threshold'ing* can be done in various ways which are further
#' described in the details.
#'
#' In case a [`SpatRaster`] is provided as input in this function for
#' \code{obj}, it is furthermore necessary to provide a [`sf`] object for
#' validation as there is no [`DistributionModel`] to read this information
#' from.
#'
#' **Note:** This of course also allows to estimate the threshold based on withheld data, for instance
#' those created from an a-priori cross-validation procedure.
#'
#' For [`BiodiversityScenario`] objects, adding this function to the processing
#' pipeline stores a threshold attribute in the created [scenario] object.
#'
#' @param obj A trained [`DistributionModel`] or alternatively a [`SpatRaster`] object.
#' @param method A specifc method for thresholding. See details for available options.
#' @param value A [`numeric`] value for thresholding if method is fixed (Default: \code{NULL}).
#' @param point A [`sf`] object containing observational data used for model training.
#' @param field_occurrence A [`character`] location of biodiversity point records.
#' @param format [`character`] indication of whether \code{"binary"}, \code{"normalize"}
#' or \code{"percentile"} formatted thresholds are to be created (Default: \code{"binary"}).
#' Also see Muscatello et al. (2021).
#' @param return_threshold Should threshold value be returned instead (Default: \code{FALSE})
#' @param ... other parameters not yet set.
#'
#' @details The following options are currently implemented:
#' * \code{'fixed'} = applies a single pre-determined threshold. Requires \code{value}
#' to be set.
#' * \code{'mtp'} = minimum training presence is used to find and set the lowest
#' predicted suitability for any occurrence point.
#' * \code{'percentile'} = For a percentile threshold. A \code{value} as parameter
#' has to be set here.
#' * \code{'min.cv'} = Threshold the raster so to minimize the coefficient of
#' variation (cv) of the posterior. Uses the lowest tercile of the cv in space.
#' Only feasible with Bayesian engines.
#' * \code{'TSS'} = Determines the optimal TSS (True Skill Statistic). Requires
#' the \code{"modEvA"} package to be installed.
#' * \code{'kappa'} = Determines the optimal kappa value (Kappa). Requires the
#' \code{"modEvA"} package to be installed.
#' * \code{'F1score'} = Determines the optimal F1score (also known as Sorensen
#' similarity). Requires the \code{"modEvA"} package to be installed.
#' * \code{'F1score'} = Determines the optimal sensitivity of presence records.
#' Requires the \code{"modEvA"} package to be installed.
#' * \code{'Sensitivity'} = Determines the optimal sensitivity of presence records.
#' Requires the \code{"modEvA"} package to be installed.
#' * \code{'Specificity'} = Determines the optimal sensitivity of presence records.
#' Requires the \code{"modEvA"} package to be installed.
#' * \code{'AUC'} = Determines the optimal AUC of presence records. Requires the
#' \code{"modEvA"} package to be installed.
#'
#' @returns A [SpatRaster] if a [SpatRaster] object as input. Otherwise the threshold
#' is added to the respective [`DistributionModel`] or [`BiodiversityScenario`] object.
#'
#' @references
#' * Lawson, C.R., Hodgson, J.A., Wilson, R.J., Richards, S.A., 2014. Prevalence,
#' thresholds and the performance of presence-absence models. Methods Ecol. Evol.
#' 5, 54–64. https://doi.org/10.1111/2041-210X.12123
#' * Liu, C., White, M., Newell, G., 2013. Selecting thresholds for the prediction
#' of species occurrence with presence-only data. J. Biogeogr. 40, 778–789. https://doi.org/10.1111/jbi.12058
#' * Muscatello, A., Elith, J., Kujala, H., 2021. How decisions about fitting
#' species distribution models affect conservation outcomes. Conserv. Biol. 35, 1309–1320.
#' https://doi.org/10.1111/cobi.13669
#'
#' @seealso \code{"modEvA"}
#'
#' @examples
#' \dontrun{
#'  # Where mod is an estimated DistributionModel
#'  tr <- threshold(mod)
#'  tr$plot_threshold()
#' }
#'
#' @name threshold
NULL

#' @rdname threshold
#' @export
methods::setGeneric(
  "threshold",
  signature = methods::signature("obj", "method", "value"),
  function(obj, method = 'mtp', value = NULL, point = NULL, field_occurrence = "observed", format = "binary", return_threshold = FALSE, ...) standardGeneric("threshold"))

#' @rdname threshold
methods::setMethod(
  "threshold",
  methods::signature(obj = "ANY"),
  function(obj, method = 'mtp', value = NULL, point = NULL, field_occurrence = "observed", format = "binary", return_threshold = FALSE, ...) {
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
                                  # modEvA measures
                                  'TSS','kappa','F1score','Sensitivity','Specificity',
                                  'Misclass','Omission','Commission','Precision',
                                  'PPI','PAI','OddsRatio'), several.ok = FALSE)

    # If method is min.cv, check that posterior is accessible
    if(method == "min.cv") assertthat::assert_that("cv" %in% names(ras),
                                                   msg = "Method min.cv requires a posterior prediction and coefficient of variation!")

    # Get all point data in distribution model
    if(is.null(point)){
      if(method %in% c('TSS','kappa','F1score','Sensitivity','Specificity',
                       'Misclass','Omission','Commission','Precision',
                       'PPI','PAI','OddsRatio')){
        if(getOption('ibis.setupmessages')) myLog('[Threshold]','yellow','Ideally thresholds are created with independent data.\n Using training data.')
      }
      point <- collect_occurrencepoints(model = model,
                                        include_absences = TRUE,
                                        point_column = field_occurrence,
                                        addName = TRUE, tosf = TRUE
                                        )
    } else {
      assertthat::assert_that(inherits(point, 'sf'),
                              utils::hasName(point, field_occurrence),
                              msg = paste0("Field ", field_occurrence, " not found in the layer!"))
      assertthat::assert_that(sf::st_crs(point) == sf::st_crs(obj$get_data('prediction')))
    }

    # If TSS or kappa is chosen, check whether there is poipa data among the sources
    if((!any(point[, field_occurrence, drop = TRUE]==0) & method %in% c('TSS','kappa','F1score','Sensitivity','Specificity'))){
      if(getOption('ibis.setupmessages')) myLog('[Threshold]','red','Threshold method needs absence-data. Generating some now...')
      bg <- terra::rasterize(obj$model$background, emptyraster(obj$get_data('prediction')))
      ass <- model$biodiversity[[1]]$pseudoabsence_settings
      if(is.null(ass)) ass <- getOption("ibis.pseudoabsence") # Get Default settings

      # Rename geometry to be consistent
      point <- rename_geometry(point, "geometry")

      suppressMessages(
        abs <- add_pseudoabsence(df = point,
                                 field_occurrence = field_occurrence,
                                 template = bg,
                                 # Assuming that settings are comparable among objects
                                 settings = ass
        )
      )
      abs <- subset(abs, select = c('x','y'));abs[, field_occurrence] <- 0
      abs <- guess_sf(abs)
      abs$name <- 'Background point'; abs$type <- "generated"
      suppressWarnings(
        abs <- sf::st_set_crs(abs, value = sf::st_crs(obj$get_data('prediction')))
      )
      point <- point |> dplyr::select(dplyr::all_of(field_occurrence), geometry, dplyr::any_of(c("name", "type")))
      abs <- abs |> dplyr::select(dplyr::all_of(field_occurrence), geometry, dplyr::any_of(c("name", "type")))
      point <- dplyr::bind_rows(point,abs)
      rm(abs)
    }

    # Convert to sf
    if(!inherits(point,"sf")){ point <- guess_sf(point) }

    # Now self call threshold
    out <- .stackthreshold(obj = ras,method = method, value = value, point = point,
                           field_occurrence = field_occurrence, format = format,
                           return_threshold = return_threshold)
    assertthat::assert_that(is.Raster(out))
    # Add result to new obj and clean up old thresholds before
    tr_lyr <- grep('threshold', obj$show_rasters(),value = TRUE)
    new_obj <- obj
    if(length(tr_lyr)>0) for(v in tr_lyr) new_obj$rm_threshold()
    new_obj <- new_obj$set_data(paste0("threshold_", method), out)
    # Return altered object
    return(new_obj)
  }
)

#' @noRd
#'
#' @keywords internal
.stackthreshold <- function(obj, method = 'fixed', value = NULL,
                            point = NULL, field_occurrence = "observed", format = "binary", return_threshold = FALSE, ...) {
  assertthat::assert_that(is.Raster(obj),
                          is.character(method),
                          inherits(point,'sf'),
                          is.character(field_occurrence),
                          is.null(value) || is.numeric(value),
                          is.character(format)
  )
  # Match format
  format <- match.arg(format, c("binary", "normalize", "percentile"), several.ok = FALSE)

  # Check for field occurrence field
  if(!is.null(point)){
    assertthat::assert_that(inherits(point, 'sf'),
                            utils::hasName(point, field_occurrence),
                            msg = paste0("Field ", field_occurrence, " not found in the layer!"))
  }

  # Apply threshold on each entry
  if(return_threshold){
    # Return the threshold directly
    out <- vector()
    for(i in names(obj)) out <- c(out, threshold(obj[[i]], method = method,
                                                                value = value, point = point, field_occurrence = field_occurrence,
                                                 format = format, return_threshold = return_threshold, ...) )
    names(out) <- names(obj)
  } else {
    # Return the raster instead
    if(method == "min.cv"){
      # If the coefficient of variation is to be minmized, mask first all values with the threshold only
      assertthat::assert_that(terra::nlyr(obj)>2, "sd" %in% names(obj))
      # Get global coefficient of variation
      errortr <- terra::global(obj[["cv"]], fun = quantile, na.rm = TRUE)[[3]]
      assertthat::assert_that(is.numeric(errortr))
      # Create mask
      mm <- obj[["cv"]]
      mm[mm > errortr] <- NA
      obj <- terra::mask(obj, mm); rm(mm)
      # Set the value to errortr
      value <- errortr
    }
    # Now loop
    out <- threshold(obj, method = method,
                     value = value, point = point, field_occurrence = field_occurrence,
                     format = format,
                     return_threshold = return_threshold, ...)
  }
  return(out)
}

#' @rdname threshold
methods::setMethod(
  "threshold",
  methods::signature(obj = "SpatRaster"),
  function(obj, method = 'fixed', value = NULL, point = NULL, field_occurrence = "observed", format = "binary", return_threshold = FALSE) {
    assertthat::assert_that(is.Raster(obj),
                            inherits(obj,'SpatRaster'),
                            is.character(method),
                            is.null(value) || is.numeric(value),
                            is.character(format),
                            is.logical(return_threshold)
    )
    # Match format
    format <- match.arg(format, c("binary", "normalize", "percentile"), several.ok = FALSE)

    # If poi is set, try to convert sf
    if(!is.null(point)) try({point <- sf::st_as_sf(point)}, silent = TRUE)
    assertthat::assert_that(is.null(point) || inherits(point,'sf'))

    # Match to correct spelling mistakes
    method <- match.arg(method, c('fixed','mtp','percentile','min.cv',
                                  # modEvA measures
                                  'TSS','kappa','F1score','Sensitivity','Specificity',
                                  'Misclass','Omission','Commission','Precision',
                                  'PPI','PAI','OddsRatio'), several.ok = FALSE)

    # Check that raster has at least a mean prediction in name
    if(!is.null(point)) {
      assertthat::assert_that(utils::hasName(point,field_occurrence),
                              msg = "Provided point data needs to include specified occurrence column!")
      # If observed is a factor, convert to numeric
      if(is.factor(point[, field_occurrence, drop = TRUE])){
        point[, field_occurrence] <- as.numeric(as.character(point[, field_occurrence, drop = TRUE]))
      }
      assertthat::assert_that(unique(sf::st_geometry_type(point)) %in% c('POINT','MULTIPOINT'))
      assertthat::assert_that(utils::hasName(point, field_occurrence))
      poi_pres <- point[point[, field_occurrence, drop = TRUE] > 0, ]  # Remove any eventual absence data for a poi_pres evaluation
    } else poi_pres <- NULL

    # Loop through each raster
    out <- c()
    for(val in names(obj)){
      # Get the raster layer
      raster_thresh <- subset(obj, val)

      # Specify by type:
      if(method == "fixed"){
        # Fixed threshold. Confirm to be set
        assertthat::assert_that(is.numeric(value), msg = 'Fixed value is missing!')
        tr <- value
      } else if(method == "mtp"){
        assertthat::assert_that(!is.null(poi_pres),msg = "Threshold method requires supplied point data!")
        # minimum training presence
        pointVals <- get_rastervalue(coords = poi_pres, env = raster_thresh)[[val]]
        # Minimum threshold
        tr <- min( stats::na.omit(pointVals) )

      } else if(method == "percentile"){
        assertthat::assert_that(!is.null(poi_pres), msg = "Threshold method requires supplied point data!")
        pointVals <- get_rastervalue(coords = poi_pres, env = raster_thresh)[[val]]
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
        assertthat::assert_that(!is.null(poi_pres),msg = "Threshold method requires supplied point data!")
        assertthat::assert_that(!is.null(value),msg = "Global minimum cv needs to be supplied as value!")
        pointVals <- get_rastervalue(coords = poi_pres, env = raster_thresh)[[val]] # Extract point only estimates

        # Get standard deviation and calculate percentile
        tr <- min( stats::na.omit(pointVals) )
        names(tr) <- "tr"
        names(value) <- "min.cv"
        # Combine as a vector
        tr <- c(tr, value)

      } else {
        # Optimized threshold statistics using the modEvA package
        # FIXME: Could think of porting these functions but too much effort for
        # now. Rather have users install the package here
        check_package("modEvA")
        # Assure that point data is correctly specified
        assertthat::assert_that(inherits(point, 'sf'), utils::hasName(point, field_occurrence))
        point[, field_occurrence] <- ifelse(point[, field_occurrence, drop = TRUE] > 1,
                                              1, point[, field_occurrence, drop = TRUE])  # Ensure that observed is <=1
        assertthat::assert_that(all( unique(point[, field_occurrence, drop = TRUE]) %in% c(0,1) ))

        # Re-extract point vals but with the full dataset
        pointVals <- get_rastervalue(coords = point, env = raster_thresh)[[val]]
        assertthat::assert_that(length(pointVals)>2)
        # Calculate the optimal thresholds
        suppressWarnings(
          opt <- modEvA::optiThresh(obs = point[, field_occurrence, drop = TRUE],
                                    pred = pointVals, measures = method,
                                    optimize = "each", plot = FALSE)
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
        o <- tr
        names(o) <- method
        out <- c(out, o)
      } else {
        # Finally threshold the raster
        # Process depending on format
        if(format == "binary"){
          # Default is to create a binary presence-absence. Otherwise truncated hinge
          raster_thresh[raster_thresh < tr[1]] <- 0
          raster_thresh[raster_thresh >= tr[1]] <- 1
          raster_thresh <- terra::as.factor(raster_thresh)
        } else if(format == "normalize"){
          raster_thresh[raster_thresh < tr[1]] <- NA
          # If truncate, ensure that resulting values are normalized
          raster_thresh <- predictor_transform(raster_thresh, option = "norm")
          raster_thresh[is.na(raster_thresh)] <- 0
          raster_thresh <- terra::mask(raster_thresh, obj[val]>=0)
          base::attr(raster_thresh, 'truncate') <- TRUE # Legacy truncate attribute

        } else if(format == "percentile") {
          raster_thresh[raster_thresh < tr[1]] <- NA
          raster_thresh <- predictor_transform(raster_thresh, option = "percentile")
          raster_thresh <- terra::mask(raster_thresh, obj[val]>=0)
          base::attr(raster_thresh, 'truncate') <- TRUE
        }

        names(raster_thresh) <- paste0('threshold_',val,'_',method)
        # Assign attributes
        base::attr(raster_thresh, 'method') <- method
        base::attr(raster_thresh, 'format') <- format
        base::attr(raster_thresh, 'threshold') <- tr

        # Append
        suppressWarnings( out <- c(out, raster_thresh) )
      }
    }

    # Return output
    if(is.list(out)) out <- do.call(c, out)
    return( out )
  }
)

#### For scenarios ####

#' Thresholds in scenario estimation
#'
#' @param obj A [BiodiversityScenario] object to which an existing threshold is
#' to be added.
#' @param value A [`numeric`] value specifying the specific threshold for scenarios
#' (Default: \code{NULL} Grab from object).
#' @param ... Any other parameter. Used to fetch value if set somehow.
#'
#' @rdname threshold
methods::setMethod(
  "threshold",
  methods::signature(obj = "BiodiversityScenario"),
  function(obj, value = NULL, ...) {

    # Get Dots
    dots <- list(...)
    # Check for value parameter in dots if tr is null
    if(is.null(value) && ("tr" %in% names(dots))){
      value <- dots[["tr"]]
      assertthat::assert_that(is.numeric(value),
                              msg = "Parameter value not found and other numeric values not found?")
    }

    # Assert that predicted raster is present
    assertthat::assert_that( is.Raster(obj$get_model()$get_data('prediction')) )

    # Unless set, check
    if(is.null(value)){
      # Check that a threshold layer is available and get the methods and data from it
      assertthat::assert_that( length( grep('threshold', obj$get_model()$show_rasters()) ) >0 ,
                               msg = 'Call \' threshold \' for prediction first!')
      # Get threshold layer
      tr_lyr <- grep('threshold', obj$get_model()$show_rasters(),value = TRUE)
      if(length(tr_lyr)>1) warning("There appear to be multiple thresholds. Using the first one.")
      ras_tr <- obj$get_model()$get_data( tr_lyr[1] )
      value <- attr(ras_tr[[1]], 'threshold')
      names(value) <- attr(ras_tr[[1]], 'method')
    } else {
      assertthat::assert_that(is.numeric(value))
      names(value) <- "fixed"
      # Check if scenario is already fitted
      proj <- obj$get_data()
      if(!is.Waiver(proj)){
        # Existing projection ?
        if(getOption('ibis.setupmessages')) myLog('[Threshold]',
                                                  'green',
                                                  'Projection found. Applying threshold!')
        new <- proj |> dplyr::select(suitability)
        obj <- obj$apply_threshold(tr = value)
        return( obj )
        invisible()
      }
    }
    bdproto(NULL, obj, threshold = value)
  }
)

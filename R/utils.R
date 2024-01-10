#' Inverse of in call for convenience
#'
#' @description
#' Calculates the set of entries not present in the second vector
#'
#' @param a First [`vector`] object.
#' @param b Second [`vector`] object.
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
`%notin%` = function(a, b){!(a %in% b)}

#' Custom messaging function for scripts
#'
#' @description This functions prints a message with a custom header and colour.
#'
#' @param title The title in the log output
#' @param col A [`character`] indicating the text colour to be used. Supported
#' are \code{'green'} / \code{'yellow'} / \code{'red'}
#' @param ... Any additional outputs or words for display
#'
#' @keywords utils
#'
#' @examples
#' myLog("[Setup]", "red", "Some error occurred during data preparation.")
#'
#' @keywords internal
myLog <- function(title = "[Processing]", col = 'green', ...) {
  assertthat::assert_that(col %in% c('green','yellow','red'))
  textwrap <- switch (col,
    'green' = text_green,
    'yellow' = text_yellow,
    'red' = text_red
  )
  message(textwrap(
    paste0(title,' ', Sys.time(), " | ", ...)
      )
    )
}

#' Colour helpers for message logs
#'
#' @param text A [`character`].
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
text_red <- function(text) { paste0('\033[31m',text,'\033[39m') }

#' @inheritParams text_red
text_yellow <- function(text) { paste0('\033[33m',text,'\033[39m') }

#' @inheritParams text_red
text_green <- function(text) { paste0('\033[32m',text,'\033[39m') }

#' Calculate the mode
#'
#' @param A [`vector`] of values or characters.
#'
#' @keywords utils
#'
#' @noRd
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Check whether function exist in name space
#'
#' @param x The [character] name of a package from which a function is needed.
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
check_package <- function(x) {
  assertthat::assert_that(is.character(x))
  if (!requireNamespace(x, quietly = TRUE)) {
    stop(paste0("Package \"",x,"\" needed for this function to work. Please install it."),
         call. = FALSE)
  }
}

#' Camel case conversion of a string
#'
#' @param x A [`vector`] or [`character`] object.
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
to_camelcase <- function(x){
  assertthat::assert_that(is.character(x) || is.vector(x))
  substr(x, 1, 1) <- toupper(
      substr(x, 1, 1)
    )
  x
}

#' Atomic representation of a name
#'
#' @description
#' Return a pretty character representation of an object with elements and names.
#'
#' @param x A [`vector`] object
#' @return [`character`] object.
#'
#' @keywords utils
#' @concept function taken from `prioritizr` package
#'
#' @examples
#' name_atomic(letters)
#' name_atomic(letters, "characters")
#'
#' @noRd
#'
#' @keywords internal
name_atomic <- function(x, description = "") {
  n <- length(x)
  if (nchar(description) > 0)
    description <- paste0(" ", description)
  if (length(x) <= 4) {
    x <- x[seq_len(min(length(x), 4))]
  } else {
    x <- c(x[seq_len(min(length(x), 3))], "...")
  }
  paste0(paste(x, collapse = ", "), " (", n, description, ")")
}

#' Aligns text with new characters
#'
#' @description
#' Format text by adding a certain number of spaces after new line characters.
#'
#' @param x [`character`] text.
#' @param n [`integer`] number of spaces.
#'
#' @return [`character`].
#'
#' @keywords utils
#' @concept function taken from `prioritizr` package
#'
#' @examples
#' # make some text
#' original_text <- "animals: horse\npig\nbear"
#'
#' # print text
#' message(original_text)
#'
#' # this look really ugly so we will align it
#' aligned_text <- align_text(original_text, 9)
#'
#' # print aligned text
#' message(aligned_text)
#'
#' @noRd
align_text <- function(x, n) {
  assertthat::assert_that(assertthat::is.string(x), assertthat::is.count(n))
  if (!grepl("\n", x))
    return(x)
  return(gsub("\n", paste0("\n", paste(rep(" ", n), collapse = "")), x,
              fixed = TRUE))
}

#' Convert character to capital text
#'
#' @param x [`character`] text.
#'
#' @keywords utils
#'
#' @examples
#' capitalize_text('presence')
#' capitalize_text('ducks are the best birds')
#'
#' @noRd
capitalize_text <- function(x) {
  assertthat::assert_that(is.character(x))
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

#' Guess time to Posix
#'
#' @description This little wrapper converts and ensures that a vector of time
#' objects are in POSIXct format.
#'
#' @param vec A [`vector`] with [`numeric`] or [`Posixct`] data
#'
#' @keywords utils
#'
#' @noRd
to_POSIXct <- function(vec){
  # Check th
  # Parse differently depending on time
  if(inherits(vec, "POSIXct")){
    out <- vec
  } else if(inherits(vec, "units") || inherits(vec, "Date")){
    check_package("units")
    # Try and format directly to posixct
    out <- as.POSIXct(vec)
    assertthat::assert_that(any(!is.na.POSIXlt(out)))
  } else if(inherits(vec, "numeric")){
    if(all(nchar(vec)==4)){
      # Assume that the numeric is a year
      vec <- paste0(vec, "-01-01")
      out <- as.POSIXct(vec)
    }
  } else if(inherits(vec, "character")){
    # Try and convert to posix directly
    out <- as.POSIXct(vec)
    if(any(is.na.POSIXlt(out))){
      # Situation not yet encountered. To be added when use cases are known.
      message("Date formats probably need some more prior handling.")
    }
  }
  return(out)
}

#' Hingeval transformation
#'
#' @param x A [`vector`] with numeric values.
#' @param min [`numeric`] minimum value for the hinge transformation.
#' @param max [`numeric`] maximum value for the hinge transformation.
#'
#' @noRd
#'
#' @keywords internal
hingeval <- function (x, min, max){
  ifelse(is.na(x),NA, pmin(1, pmax(0, (x - min)/(max - min),na.rm = TRUE),na.rm = TRUE))
}

#' Threshold transformation
#' @param x A [`vector`] with numeric values.
#' @param knot [`numeric`] threshold value as cutoff.
#'
#' @noRd
#'
#' @keywords internal
thresholdval <- function(x, knot) {
    ifelse(x >= knot, 1, 0)
}

#' Sanitize variable names
#'
#' @description Prepared covariates often have special characters in their
#' variable names which can or can not be used in formulas or cause errors for
#' certain engines. This function converts special characters of variable names
#' into a format
#'
#' @param names A [`vector`] of [`character`] vectors to be sanitized.
#'
#' @returns A [`vector`] of sanitized [`character`].
#'
#' @keywords utils
#' @concept Inspired from [`inlabru`] \code{"bru_standardise_names"} function.
#'
#' @examples
#' # Correct variable names
#' vars <- c("Climate-temperature2015", "Elevation__sealevel", "Landuse.forest..meanshare")
#' sanitize_names(vars)
#'
#' @export
sanitize_names <- function(names){
  assertthat::assert_that(
    length(names) > 0
  )
  # Convert the variable names
  new_names <- vapply(names, function(x) {
    gsub("[-() ]", "_", x = x, fixed = FALSE)
  }, "name")

  not_ok <- grepl("__", x = new_names)
  while (any(not_ok)) {
    new_names[not_ok] <- vapply(new_names[not_ok], function(x) {
      gsub("__", "_", x = x, fixed = FALSE)
    }, "name")
    not_ok <- grepl("__", x = new_names)
  }
  assertthat::assert_that(length(new_names) == length(names))
  return(
    as.character(new_names)
    )
}

#' Parallel computation of function
#'
#' @description Some computations take considerable amount of time to execute.
#' This function provides a helper wrapper for running functions of the
#' [`apply`] family to specified outputs.
#'
#' @param X A [`list`], [`data.frame`] or [`matrix`] object to be fed to a single
#' core or parallel [apply] call.
#' @param FUN A [`function`] passed on for computation.
#' @param cores A [numeric] of the number of cores to use (Default: \code{1}).
#' @param approach [`character`] for the parallelization approach taken (Options:
#' \code{"parallel"} or \code{"future"}).
#' @param export_package A [`vector`] with packages to export for use on parallel
#' nodes (Default: \code{NULL}).
#'
#' @details By default, the [parallel] package is used for parallel computation,
#' however an option exists to use the [future] package instead.
#'
#' @keywords utils
#'
#' @examples
#' \dontrun{
#'  run_par(list, mean, cores = 4)
#' }
#'
#' @noRd
run_parallel <- function(X, FUN, cores = 1, approach = "parallel", export_packages = NULL, ...) {
  assertthat::assert_that(
    is.list(X) || is.data.frame(X) || is.matrix(X),
    is.function(FUN),
    is.numeric(cores),
    is.null(export_packages) || is.character(export_packages)
  )
  # Match approach
  approach <- match.arg(approach, c("parallel", "future"), several.ok = FALSE)

  # Collect dots
  dots <- list(...)

  if(!is.list(X)){
    # Convert input object to a list of split parameters
    n_vars <- nrow(X)
    chunk_size <- ceiling(n_vars / cores)
    n_chunks <- ceiling(n_vars / chunk_size)
    chunk_list <- vector(length = n_chunks, mode = "list")

    for (i in seq_len(n_chunks)) {
      if ((chunk_size * (i - 1) + 1) <= n_vars) {
        chunk <- (chunk_size * (i - 1) + 1):(min(c(chunk_size *
                                                     i, n_vars)))
        chunk_list[[i]] <- X[chunk, ]
      }
    }
    assertthat::assert_that(sum(sapply(chunk_list, nrow)) == nrow(X))
    X <- chunk_list;rm(chunk_list)
    input_type = "data.frame" # Save to aggregate later again
  } else { input_type = "list"}

  # Process depending on cores
  if (cores == 1) {
    out <- lapply(X, FUN, ...)
  } else {
      if(approach == "parallel"){
        # check_package('doParallel')
        # require(foreach)
        # isTRUE(Sys.info()[["sysname"]] == "Windows")
        # Other operating systems
        if(!isTRUE(Sys.info()[["sysname"]] == "Windows") && is.list(X)) {
          out <- parallel::mclapply(X = X, FUN = FUN, mc.cores = cores,
                                    ...)
        } else {
          # Other operating systems
          cl <- parallel::makePSOCKcluster(cores)
          on.exit(parallel::stopCluster(cl))
          if(!is.null(export_packages)){
            # Send all specified packages to the cluster
            for(val in export_packages){
              parallel::clusterExport(cl, varlist = val,
                                      envir = as.environment(asNamespace(val)))
            }
          }
          out <- parallel::parLapply(cl = cl, X = X, fun = FUN, ...)
        }
        # out <- foreach::foreach(z = iterators::iter(X),
        #                .combine = ifelse(input_type!="list", "rbind", foreach:::defcombine),
        #                .inorder = FALSE,
        #                .multicombine = TRUE,
        #                .errorhandling = 'stop',
        #                .export = c("FUN"),
        #                .packages = export_packages,
        #                ...
        # ) %dopar% { return( FUN(z, ...) ) }
      } else {
        # Check that future is loaded
        check_package('future.apply')
        # Check that plan for future has been set up!
        assertthat::assert_that( getOption("ibis.use_future") == TRUE,
                                 msg = "Set up a future plan via [ibis_future] to use this approach.")
        out <- future.apply::future_lapply(cl = cl, X = X, fun = FUN, ...)
      }
  }
  # If input data was not a list, combine again
  if(input_type != "list" && is.list(out)){
    out <- do.call(rbind, out)
  }
  return( out )
}

#' Clamp a predictor matrix by given values
#'
#' @description To limit extreme extrapolation it is possible to \code{'clamp'}
#' an existing projection to the range of predictor values observed during model
#' training. This function takes an internal model matrix and restricts the
#' values seen in the predictor matrix to those observed during training.
#'
#' @param model A [`list`] with the input data used for inference. Created during model setup.
#' @param pred An optional [`data.frame`] of the prediction container.
#'
#' @note This function is meant to be used within a certain \code{"engine"} or
#' within [`project`].
#'
#' @returns A [`data.frame`] with the clamped predictors.
#'
#' @keywords utils
#'
#' @references
#' Phillips, S. J., Anderson, R. P., Dudík, M., Schapire, R. E., & Blair, M. E. (2017).
#' Opening the black box: An open-source release of Maxent. Ecography.
#' https://doi.org/10.1111/ecog.03049
#'
#' @noRd
#'
#' @keywords internal
clamp_predictions <- function(model, pred){
  assertthat::assert_that(
    is.list(model),
    assertthat::has_name(model, "biodiversity"),
    (is.data.frame(pred) || is.matrix(pred)) || missing(pred)
  )

  # For each biodiversity dataset, calculate the range of predictors observed
  vars_clamp <- data.frame()
  for(ds in model$biodiversity){
    # Calculate range for each variable
    rr <- apply(ds$predictors[,ds$predictors_names], 2, function(z) range(z, na.rm = TRUE)) |>
      t() |> as.data.frame() |> tibble::rownames_to_column("variable")
    names(rr) <- c("variable", "min", "max")
    vars_clamp <- rbind(vars_clamp, rr)
    rm(rr)
  }
  # Aggregate if multiple variables
  if(anyDuplicated(vars_clamp$variable)){
    o1 <- aggregate(variable ~ min, data = vars_clamp,
              FUN = function(x) min(x) )
    o2 <- aggregate(variable ~ max, data = vars_clamp,
                    FUN = function(x) max(x) )
    vars_clamp <- merge(o1,o2)
  }
  # --- #
  # Now clamp either predictors
  if(missing(pred)) pred <- model$predictors

  # Now clamp the prediction matrix with the clamped variables
  for (v in intersect(vars_clamp$variable, names(pred))) {
    pred[, v] <- pmin(
      pmax(pred[, v], vars_clamp$min[vars_clamp==v] ),
      vars_clamp$max[vars_clamp==v])
  }

  assertthat::assert_that( is.data.frame(pred) || is.matrix(pred),
                           nrow(pred)>0)
  return(pred)
}

#' Outlier detection via reverse jackknife
#'
#' @description Implementation of a Reverse Jackknife procedure as described by
#' Chapman (2005). Can be used to identify outliers in environmental predictors
#' or predictions.
#'
#' @param vals A [`numeric`] vector from which outliers are to be identified and
#' removed.
#' @param procedure [`character`] denoting what to do with outliers. Options include:
#' \code{'missing'} (Default) and \code{'remove'}, with the former replacing the
#' outliers with \code{NA} and the latter removing them.
#'
#' @keywords utils
#'
#' @references
#' * Chapman, A.D. (2005) Principles and Methods of Data Cleaning - Primary Species
#' and Species- Occurrence Data, version 1.0. Report for the Global Biodiversity
#' Information Facility, Copenhagen.
#' @source [`bioGeo`] package code served as inspiration
#'
#' @noRd
rm_outlier_revjack <- function(vals, procedure = "missing"){
  assertthat::assert_that(
    is.numeric(vals),
    length(vals)>0,
    is.character(procedure)
  )
  procedure <- match.arg(procedure, c("missing", "remove"), several.ok = FALSE)

  v2 <- vals # Make a copy
  vals <- unique(vals)
  lgh <- length(vals) - 1
  t1 <- (0.95 * sqrt(length(vals))) + 0.2
  x <- sort(vals)
  y <- rep(0, lgh)
  for (i in seq_len(lgh)) {
    x1 <- x[i + 1]
    if (x[i] < mean(vals, na.rm = TRUE)) {
      y[i] <- (x1 - x[i]) * (mean(vals, na.rm = TRUE) - x[i])
    } else {
      y[i] <- (x1 - x[i]) * (x1 - mean(vals, na.rm = TRUE))
    }
  }
  my <- mean(y, na.rm = TRUE)
  z <- y / (sqrt(sum((y - my)^2, na.rm = TRUE) / lgh))
  out <- rep(0, length(v2))
  if (any(z > t1, na.rm = TRUE)) {
    f <- which(z > t1)
    vals <- x[f]
    if (vals < stats::median(x, na.rm = TRUE)) {
      xa <- (v2 <= vals) * 1
      out <- out + xa
    }
    if (vals > stats::median(x, na.rm = TRUE)) {
      xb <- (v2 >= vals) * 1
      out <- out + xb
    }
  } else {
    out <- out
  }
  # Which ones are outliers?
  found <- which(out == 1)
  if(length(found)>0) {
    if(procedure == "missing") v2[found] <- NA else v2 <- v2[-found]
  }
  return(v2)
}

#' Aggregate count observations to a grid
#'
#' @description This function aggregates provided point data to a reference
#' grid, by, depending on the type, either counting the number of observations
#' per grid cell or aggregating them via a sum.
#'
#' @param df A [`sf`], [`data.frame`] or [`tibble`] object containing point data.
#' @param template A [`SpatRaster`] object that is aligned with the predictors.
#' @param field_occurrence A [`character`] name of the column containing the presence
#' information (Default: \code{observed}).
#'
#' @returns A [`sf`] object with the newly aggregated points.
#'
#' @noRd
#'
#' @keywords internal
aggregate_observations2grid <- function(df, template, field_occurrence = 'observed'){
  assertthat::assert_that(
    is.data.frame(df) || inherits(df, 'sf') || tibble::is_tibble(df),
    is.Raster(template),
    is.character(field_occurrence),
    assertthat::has_name(df, field_occurrence)
  )
  # Try and guess the geometry
  if(!inherits(df, 'sf')) df <- guess_sf(df)
  assertthat::assert_that(inherits(df, 'sf'), msg = "Could not convert input to sf. Prepare data first.")
  # Add coordinates if not present
  if(!assertthat::has_name(df, 'x') && !assertthat::has_name(df, 'y')) {
    df$x <- sf::st_coordinates(df[attr(df, "sf_column")])[,1]
    df$y <- sf::st_coordinates(df[attr(df, "sf_column")])[,2]
  }

  # First take presence observations and rasterize them to reduce them to a count per grid cell
  if( max(df[[field_occurrence]],na.rm = TRUE) > 1){
    # Count the sum of them
    pres <- terra::rasterize(x = df, y =  template,
                             field = field_occurrence,
                             fun = 'sum', background = 0)

  } else {
    # Simply count them
    pres <- terra::rasterize(x = df[,c("x","y")],y = template, fun = 'length', background = 0)
  }
  assertthat::assert_that(
    is.Raster(pres), is.finite( terra::global(pres, "max", na.rm = TRUE)[1,1] )
  )
  if(inherits(df, 'sf')) df <- df |> sf::st_drop_geometry()
  # Get cell ids
  ce <- terra::cellFromXY(pres, df[,c("x","y")])
  # Remove any NA if present
  if(anyNA(ce)) ce <- subset(ce, stats::complete.cases(ce))
  # Get new presence data
  obs <- cbind(
    data.frame(observed = terra::values(pres)[ce],
               terra::xyFromCell(pres, ce) # Center of cell
    )
  ) |>
    # Unique to remove any duplicate values (otherwise double counted cells)
    unique()

  # Convert to sf again
  obs <- sf::st_as_sf(obs, coords = c("x", "y"), crs = sf::st_crs(df))
  obs$x <- sf::st_coordinates(obs[attr(obs, "sf_column")])[,1]
  obs$y <- sf::st_coordinates(obs[attr(obs, "sf_column")])[,2]

  # Set CRS again
  if(is.na(sf::st_crs(obs))){
    suppressWarnings(
      obs <- sf::st_set_crs(obs, value = sf::st_crs(template))
    )
  }
  return(obs)
}

#' Get all occurrence point locations
#'
#' @description This is a small helper function that simply goes over all
#' biodiversity sets in the model object.
#' **This function is intended to only run within ibis and with the model packages created by it.**
#'
#' @param model A [`list`] object containing the biodiversity and predictor objects.
#' @param include_absences A [`logical`] of whether absences should be included
#' (Default: \code{FALSE}).
#' @param point_column [`character`] on the column with observed values.
#' @param addName [`logical`] Should the name of the feature be added
#' (Default: \code{FALSE}).
#' @param tosf [`logical`] of whether the output should be [`sf`] object
#' (Default: \code{FALSE}).
#'
#' @returns A [`matrix`] or [`sf`] object with the newly aggregated points.
#'
#' @noRd
#'
#' @keywords internal
collect_occurrencepoints <- function(model, include_absences = FALSE,
                                     point_column = "observed",
                                     addName = FALSE,
                                     tosf = FALSE){
  assertthat::assert_that(
    is.list(model),
    assertthat::has_name(model, "id"),
    assertthat::has_name(model, "biodiversity"),
    is.character(point_column),
    is.logical(addName),
    is.logical(include_absences),
    is.logical(tosf)
  )

  # Get the locations
  locs <- lapply(model$biodiversity, function(x){
                    z <- x$observations
                    # Convert to numeric just to be sure
                    if(is.factor(z[[point_column]])){
                      z[[point_column]] <- as.numeric(as.character(z[[point_column]]))
                    }

                    if(!include_absences) z <- z[z[[point_column]] > 0,]
                    if(tosf){
                      # o <- subset(
                        o <- z |> guess_sf() |>
                                  # select = c(point_column, attr(z, "sf_column"))) |>
                        rename_geometry("geometry")
                    } else {
                      o <- sf::st_coordinates(guess_sf(z)[,1:2])
                      colnames(o) <- c("x", "y")
                      o <- as.data.frame(o)
                      o[[point_column]] <- z[,point_column, drop = TRUE]
                      o[["type"]] <- x$type
                    }
                    if(addName) suppressWarnings( o$name <- x$name )
                    return(o)
                    } )
  # Combine
  locs <- do.call(rbind, locs)
  # Remove rownames
  locs <- locs |> tibble::remove_rownames()

  if(!tosf){
    assertthat::assert_that(
      is.matrix(locs) || is.data.frame(locs),
      utils::hasName(locs, point_column)
    )
  } else {
    assertthat::assert_that(inherits(locs, "sf"))
    if(is.na(sf::st_crs(locs))){
      suppressWarnings(
        locs <- locs |> sf::st_set_crs(value = sf::st_crs(model$background))
      )
    }
  }
  return(locs)
}

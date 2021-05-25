#' Inverse of in call for convenience
#' Calculates the set of entries not present in the second vector
#'
#' @param a First [`vector`] object
#' @param b Second [`vector`] object
#' @keywords internal
#' @noRd

`%notin%` = function(a, b){!(a %in% b)}

#' Custom logging function for scripts
#'
#' \code{myLog} prints a log with a custom header
#' @param title The title in the log output
#' @param ... Any additional outputs or words for display
#' @keywords internal
#' @noRd

myLog <- function(title = "[Processing]",...) {
  cat(paste0(title,' ', Sys.time(), " | ", ..., "\n"))
}

#' Colour helpers for message logs
#' @param text A [`character`]
#' @keywords internal
#' @aliases text_red
#' @noRd
text_red <- function(text) { paste0('\033[31m',text,'\033[39m') }
#' @inheritParams text_red
#' @aliases text_yellow
text_yellow <- function(text) { paste0('\033[33m',text,'\033[39m') }

#' Calculate the mode
#' @param A [`vector`] of values or characters
#' @noRd
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#' Check whether function exist in name space
#'
#' @param x The name of a package from which a function is needed
#' @examples
#'
#' @keywords internal
#' @noRd

# You need the suggested package for this function
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
#' @keywords internal
#' @noRd
to_camelcase <- function(x){
  assertthat::assert_that(is.character(x) || is.vector(x))
  substr(x, 1, 1) <- toupper(
      substr(x, 1, 1)
    )
  x
}

#' Atomic representation of a name
#'
#' Return a pretty character representation of an object with elements and
#' names.
#' Helpful function taken from `prioritizr` package
#'
#' @param x A [`vector`] object
#' @return [`character`] object.
#' @keywords internal
#' @examples
#' name_atomic(letters)
#' name_atomic(letters, "characters")
#' @noRd
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

#' Align text
#'
#' Format text by adding a certain number of spaces after new line characters.
#'
#' Helpful function taken from `prioritizr` package
#'
#' @param x [`character`] text.
#' @param n [`integer`] number of spaces.
#' @return [`character`].
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

#' Logistic (invlogit) transformation function
#' @param x A [`numeric`] value
#' @noRd
logistic <- function(a){
  #assertthat::assert_that(is.numeric(a),msg = 'Provided value not numeric')
  exp(a)/(1 + exp(a))
}
#' Logit transformation function
#' @param x A [`numeric`] value
#' @noRd
logit <- function(a){
  assertthat::assert_that(is.numeric(a),msg = 'Provided value not numeric')
  log(a/(1-a))
}

#' Convert character to formula object
#'
#' @param x [`character`] text.
#'
#' @noRd

to_formula <- function(formula){
  # Convert to formula object
  if(!is.null(formula)) {
    formula = as.formula(formula)
  } else {
    # Asign a new waiver object
    formula = new_waiver()
  }
}

#' Create formula matrix
#'
#' Function to create list of formulas with all possible combinations of variables
#' @param varnames An input [`vector`] object of variable names
#' @param response A [`character`] object describing the response
#' @param InterceptOnly Should a model with an intercept only be fitted?
#' @param special_term Any special term to add to the formulas? Default: NULL
#' @param spde_term Any spatial covariance to add to the formulas? Default: NULL
#' @param type Currently implemented are 'All' (All possible combinations) or 'forward'
#' @returns A [`list`] object with [`formula`] objects
#' @examples \dontrun{
#' formula_combinations(form)
#' }
#' @noRd

formula_combinations <- function(varnames, response = 'Observed', InterceptOnly = TRUE,
                                 special_term = NULL, spde_term = NULL,
                                 type= 'forward'){
  assertthat::assert_that(is.vector(varnames),
                          is.character(response),
                          'purrr' %in% loadedNamespaces())
  if(!is.null(special_term)) varnames <- c(varnames, special_term)

  # Formula length
  fl <- length(varnames)

  if(tolower(type) == 'forward'){
    form_temp <- c()
    for(i in 1:fl) {
      new <- paste0(response, '~ 0 + intercept +', paste(varnames[1:i],collapse = ' + ') )
      if(!is.null(spde_term)) new <- paste0(new, ' + ',spde_term)
      form_temp <- c(form_temp, new)
      }

  } else if(tolower(type) == 'all'){
    # Construct unique combinations
    varnames_comb <- 1:length(varnames) %>%
      purrr::map(~ combn(varnames, .x) %>% apply(2, list) %>% unlist(recursive = F)) %>%
      unlist(recursive = F)

    form_temp <- varnames_comb %>% purrr::map(~paste0(response, " ~ ", paste(.x, collapse = " + ")) %>% as.formula)
  }

  if(InterceptOnly) form_temp <- form_temp %>% append(paste0(response, " ~ 1") %>% as.formula %>% list, .)

  return(form_temp)
}

#' Filter a set of correlated predictors to fewer ones
#'
#' Code inspired from the `caret` package
#'
#' @param env A [`data.frame`] with extracted environmental covariates for a given species
#' @param keep A [`vector`] with variables to keep regardless
#' @param cutoff A [`numeric`] variable specifying the maximal correlation cutoff
#' @param method Which method to use for constructing the correlation matrix (pearson|spearman|kendal)
#' @returns vector of variable names to exclude

find_correlated_predictors <- function( env, keep = NULL, cutoff = 0.9, method = 'pearson' ){
  # Security checks
  assertthat::assert_that(is.data.frame(env),
                          is.character(method),
                          is.numeric(cutoff),
                          is.null(keep) || is.vector(keep)
  )
  if(!is.null(keep)) x <- env %>% dplyr::select(-keep) else x <- env

  # Removing non-numeric columns
  non.numeric.columns <- colnames(x)[!sapply(x, is.numeric)]
  x <- x[, !(colnames(x) %in% non.numeric.columns)]

  # Get all variables that are singular or unique in value
  singular_var <- which(round( apply(x, 2, var),4) == 0)
  if(length(singular_var)>0) x <- x[,-singular_var]

  # Calculate correlation matrix
  cm <- cor(x, method = method)

  # Copied from the \code{caret} package to avoid further dependencies
  if (any(!complete.cases(cm))) stop("The correlation matrix has some missing values.")
  averageCorr <- colMeans(abs(cm))
  averageCorr <- as.numeric(as.factor(averageCorr))
  cm[lower.tri(cm, diag = TRUE)] <- NA

  # Determine combinations over cutoff
  combsAboveCutoff <- which(abs(cm) > cutoff)
  colsToCheck <- ceiling(combsAboveCutoff/nrow(cm))
  rowsToCheck <- combsAboveCutoff%%nrow(cm)

  # Exclude columns with variables over average correlation
  colsToDiscard <- averageCorr[colsToCheck] > averageCorr[rowsToCheck]
  rowsToDiscard <- !colsToDiscard

  # Get columns to discard
  deletecol <- c(colsToCheck[colsToDiscard], rowsToCheck[rowsToDiscard])
  deletecol <- unique(deletecol)

  # Which variables to discard
  o <- names(env)[deletecol]
  if(length(singular_var)>0) o <- unique( c(o,  names(singular_var) ) )
  o
}

#' Create pseudo absence points over a raster dataset
#'
#' @param env An environmental dataset either in [`data.frame`] or [`raster`] format
#' @param presence Presence records. Necessary to avoid sampling pseudo-absences over existing presence records
#' @param template If template is not null then env needs to be a [`raster`] dataset
#' @param npoints A [`numeric`] number of pseudo-absence points to create
#' @param replace Sample with replacement? (Default: False)
#' @returns A [`data.frame`] containing the newly created pseudo absence points
# TODO: Allow option to supply a bias raster for weighted sampling

create_pseudoabsence <- function(env, presence, template = NULL, npoints = 1000, replace = FALSE){
  assertthat::assert_that(
    inherits(env,'Raster') || inherits(env, 'data.frame') || inherits(env, 'tibble'),
    is.data.frame(presence),
    hasName(presence,'x') && hasName(presence,'y'),
    is.numeric(npoints),
    is.null(template) || inherits(template,'Raster')
  )
  if(is.null(template)) {
    assertthat::assert_that(inherits(env,'Raster'), msg = 'Supply a template raster or a raster file')
  }

  if(inherits(env,'Raster')) env <- as.data.frame(env, xy = TRUE)

  # Rasterize the presence estimates
  bg1 <- raster::rasterize(presence[,c('x','y')], template, fun = 'count', background = 0)
  bg1 <- raster::mask(bg1, template)

  assertthat::assert_that(
    is.finite(raster::cellStats(bg1,'max',na.rm = T))
  )
  # Generate pseudo absence data
  # Now sample from all cells not occupied
  abs <- sample(which(bg1[]==0), size = npoints, replace = replace)
  # Now get absence environmental data
  abs <- get_ngbvalue(
    coords = raster::xyFromCell(bg1, abs),
    env = env,
    field_space = c('x','y'),
    longlat = raster::isLonLat(template)
  )

  # Remove NA data in case points were sampled over non-valid regions
  abs <- subset(abs, complete.cases(abs))
  assertthat::assert_that( all( names(abs) %in% names(env) ) )

  return(abs)
}

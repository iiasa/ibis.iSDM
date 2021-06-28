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
#' @param form An input [`formula`] object
#' @param response A [`character`] object giving the response. (Default: NULL)
#' @param type Currently implemented are 'inla' (variable groups), 'All' (All possible combinations) or 'forward'
#' @returns A [`vector`] object with [`formula`] objects
#' @examples \dontrun{
#' formula_combinations(form)
#' }
#' @noRd

formula_combinations <- function(form, response = NULL, type= 'forward'){
  assertthat::assert_that(is.formula(form),
                          is.character(response) || is.null(response),
                          tolower(type) %in% c('inla','forward','all'))
  # --- #
  # Response
  if(is.null(response)) response <- all.vars(form)[1]
  # Formula terms
  te <- attr(stats::terms.formula(form),'term.label')
  # Varnames
  varnames <- all.vars(form)
  varnames <- varnames[varnames %notin% c('spde','spatial.field','observed','intercept')] # Exclude things not necessarily needed in there
  # Variable length
  fl <- length(varnames)
  # --- #
  assertthat::assert_that(fl>0, !is.null(response))

  if(tolower(type) == 'inla'){
    # INLA modelling groups
    # Instead of selecting variables piece by piece, consider individual groups
    form_temp <- c()
    val_int <- grep(pattern = 'intercept',x = te, value = T)
    val_lin <- grep(pattern = 'linear',x = te, value = T)
    val_rw1 <- grep(pattern = 'rw1',x = te,value = TRUE)
    # Alternative quadratic variables in case rw1 fails
    if(length(val_rw1)>0){
      val_quad <- all.vars(as.formula(paste('observed ~ ', paste0(val_rw1,collapse = '+'))))[-1]
    } else { val_quad <- all.vars(as.formula(paste('observed ~ ', paste0(val_lin,collapse = '+'))))[-1] }
    val_spde <- grep(pattern = 'spde',x = te,value = TRUE)
    val_ofs <- grep(pattern = 'offset',x = te,value = TRUE)

    # Construct formulas ---
    # Original form
    form_temp <- c(form_temp, deparse1(form))

    # Intercept only
    form_temp <- c(form_temp,
                   paste0(response,' ~ 0 +', paste(val_int,collapse = ' + ') ))

    # Add all linear variables as base
    form_temp <- c(form_temp,
                   paste0(response,' ~ 0 +', paste(val_int,collapse = ' + '),
                          '+', paste0(varnames, collapse = ' + ') ))

    # Intercept + linear effect
    if(length(val_lin)>0){
      form_temp <- c(form_temp,
                     paste0(response,' ~ 0 + ', paste(val_int,collapse = ' + '),
                     '+',
                     paste(val_lin,collapse = ' + '))
      )
    }
    # Intercept + rw1 effects (if existing)
    if(length(val_rw1)>0){
      form_temp <- c(form_temp,
                     paste0(response,' ~ 0 + ', paste(val_int,collapse = ' + '),
                     '+',
                     paste(val_rw1,collapse = ' + '))
      )
    }
    # Alternative formulation using quadratic
    form_temp <- c(form_temp,
                   paste0(response,' ~ 0 + ', paste(val_int,collapse = ' + '),
                          '+',
                          paste0('I(',val_quad,'^2)',collapse = ' + '))
    )

    # Intercept + spde
    if(length(val_spde)>0){
      form_temp <- c(form_temp,
                     paste0(response,' ~ 0 + ', paste(val_int,collapse = ' + '),
                     '+',
                     paste(val_spde,collapse = ' + '))
      )
    }
    # Intercept + linear + spde
    if(length(val_spde)>0 && length(val_lin)>0){
      form_temp <- c(form_temp,
                     paste0(response,' ~ 0 + ', paste(val_int,collapse = ' + '),
                     '+',
                     paste(val_spde,collapse = ' + '),'+',paste(val_lin,collapse = ' + '))
      )
      form_temp <- c(form_temp,
                     paste0(response,' ~ 0 +', paste(val_int,collapse = ' + '),
                            '+', paste0(varnames, collapse = ' + '),
                            '+',paste(val_spde,collapse = ' + ')))

    }
    # intercept + rw1 + spde
    if(length(val_spde)>0 && length(val_rw1)>0){
      form_temp <- c(form_temp,
                     paste0(response,' ~ 0 + ', paste(val_int,collapse = ' + '),
                     '+',
                     paste(val_spde,collapse = ' + '),'+',paste(val_rw1,collapse = ' + '))
      )
    }
    if(length(val_spde)>0){
      # Quad replacement
      form_temp <- c(form_temp,
                     paste0(response,' ~ 0 + ', paste(val_int,collapse = ' + '),
                            '+',
                            paste(val_spde,collapse = ' + '),'+',paste0('I(',val_quad,'^2)',collapse = ' + '))
      )
    }
    # intercept + linear + rw1 + spde
    if(length(val_rw1)>0 && length(val_lin)>0 && length(val_spde)>0){
      form_temp <- c(form_temp,
                     paste0(response,' ~ 0 + ', paste(val_int,collapse = ' + '),
                            '+',
                            paste(val_lin,collapse = ' + '),'+',paste(val_rw1,collapse = ' + '),'+',paste(val_spde,collapse = ' + '))
      )

    }
    if(length(val_spde)>0){
      # Quad replacement
      form_temp <- c(form_temp,
                     paste0(response,' ~ 0 + ', paste(val_int,collapse = ' + '),
                            '+',
                            paste(val_lin,collapse = ' + '),'+',paste0('I(',val_quad,'^2)',collapse = ' + '),'+',paste(val_spde,collapse = ' + '))
      )
    }
    # intercept + linear + offset
    if(length(val_lin)>0 && length(val_ofs)>0){
      form_temp <- c(form_temp,
                     paste0(response,' ~ 0 + ', paste(val_int,collapse = ' + '),
                            '+',
                            paste(val_lin,collapse = ' + '),'+',paste(val_ofs,collapse = ' + '))
      )
    }
    # intercept + linear + rw1 + offset
    if(length(val_rw1)>0 && length(val_lin)>0 && length(val_ofs)>0){
      form_temp <- c(form_temp,
                     paste0(response,' ~ 0 + ', paste(val_int,collapse = ' + '),
                            '+',
                            paste(val_lin,collapse = ' + '),'+',paste(val_rw1,collapse = ' + '),'+',paste(val_ofs,collapse = ' + '))
      )
    }
    if(length(val_lin)>0 && length(val_ofs)>0){
    # Quad replacement
    form_temp <- c(form_temp,
                   paste0(response,' ~ 0 + ', paste(val_int,collapse = ' + '),
                          '+',
                          paste(val_lin,collapse = ' + '),'+',
                          paste0('I(',val_quad,'^2)',collapse = ' + '),'+',paste(val_ofs,collapse = ' + '))
      )
    }

  # Other types of variable selection
  } else if(tolower(type) == 'forward'){
    # Forward variable addition
    # Note this ignores unique combinations
    form_temp <- c()
    for(i in 1:fl) {
      new <- paste0(response, '~ 0 + ',paste(val_int,collapse = '+'),'+',
                    paste(varnames[1:i],collapse = ' + ') )
      form_temp <- c(form_temp, new)
    }

  } else if(tolower(type) == 'all'){
    assertthat::assert_that('purrr' %in% loadedNamespaces())
    # Construct all possible unique combinations
    varnames_comb <- 1:length(varnames) %>%
      purrr::map(~ combn(varnames, .x) %>% apply(2, list) %>% unlist(recursive = F)) %>%
      unlist(recursive = F)

    form_temp <- varnames_comb %>% purrr::map(~paste0(response, " ~ ",
                                                      paste(val_int,collapse = '+'),'+',
                                                      paste(.x, collapse = " + ")) )
  }

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

find_correlated_predictors <- function( env, keep = NULL, cutoff = 0.7, method = 'pearson'){
  # Security checks
  assertthat::assert_that(is.data.frame(env),
                          is.character(method),
                          is.numeric(cutoff),
                          is.null(keep) || is.vector(keep)
  )
  keep <- keep[keep %in% names(env)] # Remove those not in the data.frame. For instance if a spatial effect is selected
  if(!is.null(keep) || length(keep) == 0) x <- env %>% dplyr::select(-keep) else x <- env

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

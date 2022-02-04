#' Inverse of in call for convenience
#' Calculates the set of entries not present in the second vector
#'
#' @param a First [`vector`] object.
#' @param b Second [`vector`] object.
#' @keywords internal, utils
#' @noRd
`%notin%` = function(a, b){!(a %in% b)}

#' Custom messaging function for scripts
#'
#' @description
#' This functions prints a message with a custom header and colour.
#' @param title The title in the log output
#' @param col A [`character`] indicating the text colour to be used. Supported are 'green' / 'yellow' / 'red'
#' @param ... Any additional outputs or words for display
#' @examples
#' myLog("[Setup]", "red", "Some error occurred during data preparation.")
#' @keywords internal, utils
#' @export
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
#' @param text A [`character`].
#' @keywords internal, utils
#' @aliases text_red
#' @noRd
text_red <- function(text) { paste0('\033[31m',text,'\033[39m') }
#' @inheritParams text_red
#' @aliases text_yellow
text_yellow <- function(text) { paste0('\033[33m',text,'\033[39m') }
#' @inheritParams text_red
#' @aliases text_green
text_green <- function(text) { paste0('\033[32m',text,'\033[39m') }

#' Calculate the mode
#' @param A [`vector`] of values or characters.
#' @keywords utils
#' @noRd
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#' Check whether function exist in name space
#'
#' @param x The [character] name of a package from which a function is needed.
#' @keywords internal, utils
#' @noRd
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
#' @keywords internal, utils
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
#' @param x A [`vector`] object
#' @return [`character`] object.
#' @concept function taken from `prioritizr` package
#' @keywords internal, utils
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

#' Aligns text with new characters
#'
#' Format text by adding a certain number of spaces after new line characters.
#'
#' @param x [`character`] text.
#' @param n [`integer`] number of spaces.
#' @return [`character`].
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
#' @keywords utils
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
#' @examples
#' capitalize_text('presence')
#' capitalize_text('ducks are the best birds')
#'
#' @keywords utils
#' @noRd
capitalize_text <- function(x) {
  assertthat::assert_that(is.character(x))
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

#' Convert character to formula object
#'
#' @param x [`character`] text.
#' @keywords utils
#' @noRd
to_formula <- function(formula){
  # Convert to formula object
  if(!is.null(formula)) {
    formula = as.formula(formula)
  } else {
    # Asign a new waiver object
    formula = new_waiver()
  }
  return(formula)
}

#' Parallel computation of function
#'
#' @description
#' Some computations take considerable amount of time to execute. This
#' function provides a helper wrapper for running functions of the [`apply`]
#' family to specified outputs.
#' @details
#' By default, the [parallel] package is used for parallel computation,
#' however an option exists to use the [future] package instead.
#' @param X A [`list`], [`data.frame`] or [`matrix`] object to be fed to a single core or parallel [apply] call.
#' @param FUN A [`function`] passed on for computation.
#' @param cores A [numeric] of the number of cores to use (Default: \code{1}).
#' @param approach [`character`] for the parallelization approach taken (Options: \code{"parallel"} or \code{"future"}).
#' @param export_package A [`vector`] with packages to export for use on parallel nodes (Default: \code{NULL}).
#' @examples
#' \dontrun{
#'  run_par(list, mean, cores = 4)
#' }
#' @keywords utils
#' @noRd
run_parallel <- function (X, FUN, cores = 1, approach = "parallel", export_packages = NULL, ...) {
  assertthat::assert_that(
    is.list(X) || is.data.frame(X) || is.matrix(X),
    is.function(FUN),
    is.numeric(cores),
    is.null(export_packages) || is.character(export_packages)
  )
  # Match approach
  approach <- match.arg(approach, c("approach", "future"), several.ok = FALSE)

  # Collect dots
  dots <- list(...)

  if(!is.list(X)){
    # Convert input object to a list of split parameters
    n_vars <- length( ncol(X) )
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
    X <- chunk_list
    input_type = "data.frame" # Save to aggregate later again
  } else { input_type = "list"}

  # Determine the type of input and which function to use
  if(approach == "parallel"){
    # Use parallel package
    check_package('parallel')
    parfun <- switch (class(X),
      "list" = parallel::parLapply,
      "data.frame" = parallel::parApply,
      "matrix" = parallel::parApply
    )
    if(class(X) %in% c("data.frame", "matrix")) dots[["MARGIN"]] <- 1
    if (isTRUE(Sys.info()[["sysname"]] == "Windows") && !is.list(X) ){
      # Use future instead for windows
      return(
        run_parallel(X, FUN, cores = 1, approach = "future", ...)
      )
    }
  } else {
    # Use future package
    check_package('future.apply')
    # Check that plan for future has been set up!
    assertthat::assert_that( getOption("ibis.use_future") == TRUE,
                             msg = "Set up a future plan via [ibis_future] to use this approach.")
    parfun <- switch (class(X),
                      "list" = future.apply::future_lapply,
                      "data.frame" = future.apply::future_apply,
                      "matrix" = future.apply::future_apply
    )
    if(class(X) %in% c("data.frame", "matrix")) dots[["MARGIN"]] <- 1
  }
  assertthat::assert_that(is.function(parfun))

  # Process depending on cores
  if (cores == 1) {
    out <- parfun(X, FUN, ...)
  } else {
      if(approach == "parallel"){
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
              parallel::clusterExport(cl, varlist = package_function_names(val),
                                      envir = as.environment(asNamespace(val)))
            }
          }
          out <- parfun(cl = cl, X = X, fun = FUN, ...)
        }
      } else {
        # Check that future is loaded
        out <- parfun(cl = cl, X = X, fun = FUN, ...)
      }
  }
  # If input data was not a list, combine again
  if(input_type != "list"){
    out <- do.call(rbind, out)
  }
  return( out )
}
#' Create formula matrix
#'
#' Function to create list of formulas with all possible combinations of variables
#' @param form An input [`formula`] object.
#' @param response A [`character`] object giving the response. (Default: \code{NULL})
#' @param type Currently implemented are \code{'inla'} (variable groups),
#'  \code{'All'} (All possible combinations) or \code{'forward'}.
#' @returns A [`vector`] object with [`formula`] objects.
#' @examples \dontrun{
#' formula_combinations(form)
#' }
#' @keywords utils
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
#' @param env A [`data.frame`] with extracted environmental covariates for a given species.
#' @param keep A [`vector`] with variables to keep regardless.
#' @param cutoff A [`numeric`] variable specifying the maximal correlation cutoff.
#' @param method Which method to use for constructing the correlation matrix (Options: \code{'pearson'}| \code{'spearman'}| \code{'kendal'})
#' @concept Code inspired from the [`caret`] package
#' @keywords utils
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

#' Apply the adaptive best subset selection framework on a set of predictors
#'
#' @description
#' This is a wrapper function to fit the adaptive subset selection procedure outlined
#' in Zhu et al. (2021) and Zhu et al. (2020).
#' @param env A [`data.frame`] with extracted environmental covariates for a given species.
#' @param observed A [`vector`] with the observed response variable.
#' @param family A [`character`] indicating the family the observational data originates from.
#' @param tune.type [`character`] indicating the type used for subset evaluation.
#' Options are \code{c("gic", "ebic", "bic", "aic", "cv")} as listed in [abess].
#' @param lambda A [`numeric`] single lambda value for regularized best subset selection (Default: \code{0}).
#' @param weight Observation weights. When weight = \code{NULL}, we set weight = \code{1} for each observation as default.
#' @param keep A [`vector`] with variables to keep regardless (Default: \code{NULL}).
#' @references
#' * abess: A Fast Best Subset Selection Library in Python and R. Jin Zhu, Liyuan Hu, Junhao Huang, Kangkang Jiang, Yanhang Zhang, Shiyun Lin, Junxian Zhu, Xueqin Wang (2021). arXiv preprint arXiv:2110.09697.
#' * A polynomial algorithm for best-subset selection problem. Junxian Zhu, Canhong Wen, Jin Zhu, Heping Zhang, Xueqin Wang. Proceedings of the National Academy of Sciences Dec 2020, 117 (52) 33117-33123; doi: 10.1073/pnas.2014241117
#' @keywords utils
#' @returns vector of variable names to exclude
find_subset_of_predictors <- function( env, observed, family, tune.type = "cv", lambda = 0,
                                       weight = NULL, keep = NULL){
  # Security checks
  assertthat::assert_that(is.data.frame(env),
                          is.vector(observed),
                          is.numeric(lambda),
                          is.character(tune.type),
                          is.null(weight) || is.vector(weight)
  )
  assertthat::assert_that(
    length(observed) == nrow(env), msg = "Number of observation unequal to number of covariate rows."
  )
  # Match family and type
  family <- match.arg(family, c("gaussian", "binomial", "poisson", "cox", "mgaussian", "multinomial",
                               "gamma"), several.ok = FALSE)
  tune.type <- match.arg(tune.type, c("gic", "ebic", "bic", "aic", "cv"), several.ok = FALSE)

  # Check that abess package is available
  check_package("abess")
  if(!isNamespaceLoaded("abess")) { attachNamespace("abess");requireNamespace('abess') }

  # Build model
  abess_fit <- abess::abess(x = env,
                            y = observed,
                            family = family,
                            tune.type = tune.type,
                            weight = weight,
                            always.include = keep,
                            nfolds = 100, # Increase from default 5
                            num.threads = 0
                          )

  if(anyNA(coef(abess_fit)[,1]) ) {
    # Refit with minimum support size
    abess_fit <- abess::abess(x = env,
                              y = observed,
                              family = family,
                              tune.type = tune.type,
                              weight = weight,
                              always.include = keep,
                              nfolds = 100, # Increase from default 5
                              # Minimum support site of 10% of number of covariates
                              support.size = ceiling(ncol(env) * 0.1),
                              num.threads = 0
    )

  }
  # Get best vars
  co <- coef(abess_fit, support.size = abess_fit[["best.size"]])
  co <- names( which(co[,1] != 0))
  co <- co[grep("intercept", co, ignore.case = TRUE, invert = TRUE)]
  # Make some checks on the list of reduced variables
  if(length(co) <= 2) {
    warning("Abess was likely to rigours. Likely to low signal-to-noise ratio.")
    return(NULL)
  } else {
    co
  }
}


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
              parallel::clusterExport(cl, varlist = package_function_names(val),
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
  varnames <- varnames[varnames %notin% c('spde','spatial.field','observed','Intercept')] # Exclude things not necessarily needed in there
  # Variable length
  fl <- length(varnames)
  # --- #
  assertthat::assert_that(fl>0, !is.null(response))

  if(tolower(type) == 'inla'){
    # INLA modelling groups
    # Instead of selecting variables piece by piece, consider individual groups
    form_temp <- c()
    val_int <- grep(pattern = 'Intercept',x = te, value = T)
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

#' Outlier detection via reverse jackknife
#'
#' @description
#' Implemententation of a Reverse Jackknife procedure as described by Chapman (2005).
#' Can be used to identify outliers in environmental predictors or predictions.
#' @param vals A [`numeric`] vector from which outliers are to be identified and removed.
#' @param procedure [`character`] denoting what to do with outliers.
#' Options include: \code{'missing'} (Default) and \code{'remove'}, with the former replacing the outliers with \code{NA} and the latter removing them.
#' @references
#' * Chapman, A.D. (2005) Principles and Methods of Data Cleaning - Primary Species and Species- Occurrence Data, version 1.0. Report for the Global Biodiversity Information Facility, Copenhagen.
#' @source [`bioGeo`] package code served as inspiration
#' @keywords utils
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
    if (vals < median(x, na.rm = TRUE)) {
      xa <- (v2 <= vals) * 1
      out <- out + xa
    }
    if (vals > median(x, na.rm = TRUE)) {
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
  co <- co[grep("Intercept", co, ignore.case = TRUE, invert = TRUE)]
  # Make some checks on the list of reduced variables
  if(length(co) <= 2) {
    warning("Abess was likely to rigours. Likely to low signal-to-noise ratio.")
    return(NULL)
  } else {
    co
  }
}

#' Specify settings for pseudo-absence sampling from a background
#'
#' @description
#' This function defines the settings for pseudo-absence background sampling. For many
#' engines background points are necessary to model poisson (or binomial) distributed data.
#' @details
#' Mandatory and possible parameters are:
#'
#' * \code{background} Specifies the extent over which background points are to be sampled.
#' * \code{nrpoints} The number of background points to be used.
#' * \code{method} The specific method on how pseudo-absence points should be generated. Available options
#' are \code{'random'} (Default), \code{'buffer'} to generate only within a buffer of existing points, or \code{'mcp'} to only generate within the
#' for a minimum convex polygon of provided points so that absence points are only sampled outside of it.
#' * \code{buffer_distance} The numeric value indicating the distance from presence points where absence points are to be created.
#' * \code{min_ratio} The minimum ratio of absence points relative presence points. Ensures a minimum number of
#' absence points.
#' * \code{bias} An optional bias layer over which points should preferentially be sampled.
#' @param background A [`RasterLayer`] or [`sf`] object over which background points can be sampled. Default is
#' \code{NULL} is the default and the background is then added when the sampling is first called.
#' @param nrpoints A [`numeric`] given the number of background points to be created (Default: \code{10 000}).
#' @param method [`character`] denoting how the sampling should be done. For details for options (Default: \code{"random"}).
#' @param buffer_distance [`numeric`] A distance from the observations in which pseudo-absence points are not to be generated.
#' The units in \code{m}.
#' @param mcp_inside A [`logical`] value of whether absence points should be sampled inside (Default) or outside the
#' minimum convex polygon when this is chosen (parameter \code{method = "mcp"}).
#' @param min_ratio A [`numeric`] with the minimum Ratio of background points relative to presence points.
#' Usually ignored unless the ratio exceeds the \code{nrpoints} parameters (Default: \code{0.25}).
#' @param bias A [`RasterLayer`] with the same extent and projection and background. Background points will
#' be preferentially sampled in areas with higher (!) bias. (Default: \code{NULL}).
#' @param ... Any other settings to be added to the pseudoabs settings.
#' @examples
#' \dontrun{
#' # It is also possible to match the number of presence-absence points directly.
#' pseudoabs_settings(nrpoints = 0, min_ratio = 1)
#' }
#' @name pseudoabs_settings
#' @aliases pseudoabs_settings
#' @keywords train
#' @exportMethod pseudoabs_settings
#' @export
NULL
methods::setGeneric("pseudoabs_settings",
                    signature = methods::signature("background"),
                    function(background = NULL, nrpoints = 10000, min_ratio = 0.25,
                             method = "random", buffer_distance = 10000, mcp_inside = TRUE,
                             bias = NULL, ...) standardGeneric("pseudoabs_settings"))

#' @name pseudoabs_settings
#' @rdname pseudoabs_settings
#' @usage \S4method{pseudoabs_settings}{ANY, numeric, numeric, character, numeric, logical, ANY}(background, nrpoints, min_ratio, method, buffer_distance, mcp_inside, bias)
methods::setMethod(
  "pseudoabs_settings",
  methods::signature(background = "ANY"),
  function(background = NULL, nrpoints = 10000, min_ratio = 0.25,
           method = "random", buffer_distance = 10000, mcp_inside = TRUE,
           bias = NULL, ...){
    # Check inputs
    assertthat::assert_that(
      is.Raster(background) || inherits(background, 'sf') || is.null(background),
      is.numeric(nrpoints),
      is.numeric(min_ratio),
      is.character(method),
      is.numeric(buffer_distance),
      is.Raster(bias) || is.null(bias)
    )
    method <- match.arg(method, c("random", "buffer", "mcp"), several.ok = FALSE)
    # Create the settings object
    settings <- bdproto(NULL, Settings)
    settings$name <- "Background"
    settings$set('background', background)
    # Set all options
    settings$set('nrpoints', nrpoints)
    settings$set('min_ratio', min_ratio)
    settings$set('method', method)
    settings$set('buffer_distance', buffer_distance)
    settings$set('bias', bias)
    # Other settings
    mc <- match.call(expand.dots = FALSE)
    settings$data <- c( settings$data, mc$... )

    return(settings)
  }
)

#' Add pseudo-absence points to a point data set
#'
#' @description
#' For most engines, background or pseudo-absence points are necessary. The distinction
#' lies in how the absence data are handled. For [`poisson`] distributed responses,
#' absence points are considered background points over which the intensity of sampling lambda
#' is integrated (in a classical Poisson-Process-Model).
#' In contrast in [`binomial`] distributed responses, the absence information is assumed to
#' be an adequate representation of the true absences and treated by the model as such..
#' @details
#' A [`pseudoabs_settings()`] object can be added to setup specific absence-sampling. A \code{bias} parameter
#' can be set to specify a bias layer to sample from, for instance a layer of accessibility.
#' Note that when modelling several datasets, it might make sense to check across all datasets
#' whether certain areas are truly absent.
#' By default, the pseudo-absence points are not sampled in areas in which there are already presence points.
#' @note
#' This method removes all columns from the input \code{df} object other than the \code{field_occurrence} column
#' and the coordinate columns (which will be created if not already present).
#' @param df A [`sf`], [`data.frame`] or [`tibble`] object containing point data.
#' @param field_occurrence A [`character`] name of the coloumn containing the presence information (Default: \code{observed}).
#' @param template A [`RasterLayer`] object that is aligned with the predictors (Default: \code{NULL}). If set to \code{NULL},
#' then \code{background} in the [`pseudoabs_settings()`] has to be a [`RasterLayer`] object.
#' @param settings A [`pseudoabs_settings()`] objects. Takes the default ibis settings if not set.
#' @references
#' * Stolar, J., & Nielsen, S. E. (2015). Accounting for spatially biased sampling effort in presence‐only species distribution modelling. Diversity and Distributions, 21(5), 595-608.
#' @keywords train
#' @returns A [`data.frame`] containing the newly created pseudo absence points.
#' @export
add_pseudoabsence <- function(df, field_occurrence = "observed", template = NULL, settings = getOption("ibis.pseudoabsence")){
  assertthat::assert_that(
    is.data.frame(df) || inherits(df, 'sf') || tibble::is_tibble(df),
    is.Raster(template) || is.null(template),
    is.character(field_occurrence),
    assertthat::has_name(df, field_occurrence),
    inherits(settings, "Settings")
  )
  # Check that no 0 are present, otherwise raise a warning.
  assertthat::see_if( any(df[[field_occurrence]] != 0) )

  # Try and guess the geometry
  if(!inherits(df, 'sf')) df <- guess_sf(df)
  assertthat::assert_that(inherits(df, 'sf'), msg = "Could not convert input to sf. Prepare data first.")
  # Add coordinates if not present
  if(!assertthat::has_name(df, 'x') && !assertthat::has_name(df, 'y')) {
    df$x <- sf::st_coordinates(df[attr(df, "sf_column")])[,1]
    df$y <- sf::st_coordinates(df[attr(df, "sf_column")])[,2]
  }
  # Select relevant columns and assign type
  df$type <- "Presence"

  # Check whether the background is set and if not, use the bbox
  if(is.Waiver(settings$get("background"))){
    # Check whether temlate wasn't provided
    if(!is.null(template)){
      background <- template
    } else {
      background <- sf::st_as_sf(
        sf::st_as_sfc(sf::st_bbox(df))
      );background$bg <- 1
    }
  } else {
    background <- settings$get("background")
  }
  # Check that background is a raster, otherwise rasterize with identical resolution
  if(!is.Raster(background)){
    assertthat::assert_that(is.Raster(template),
                            msg = "No suitable RasterLayer was provided through Settings or as template!")
    if("fasterize" %in% installed.packages()[,1]){
      background <- fasterize::fasterize(sf = background, raster = emptyraster(template), field = NULL)
    } else {
      background <- raster::rasterize(background, emptyraster(template), field = 1)
    }
    assertthat::assert_that(is.Raster(background))
  }

  # --- #
  # Now depending on the settings create absence points
  # Get number of points to sample and ratio
  nrpoints <- settings$get('nrpoints')
  min_ratio <- settings$get('min_ratio')

  method <- settings$get('method')
  buffer_distance <- settings$get('buffer_distance')

  # If the nr of points is 0, set it equal to the number of min_ratio or presented presence points
  if(nrpoints == 0) nrpoints <- round( nrow(df) * min_ratio )
  if(nrpoints > raster::ncell(background)) nrpoints <- raster::ncell(background)

  if(!is.Waiver(settings$get("bias"))){
    bias <- settings$get("bias")
    if(!compareRaster(bias, background, stopiffalse = FALSE)){
      #raster::compareCRS(bias, background) # Raise error if projection is different
      # Resample to ensure same coverage
      bias <- raster::crop(bias, raster::extent(background))
      bias <- raster::resample(bias, background, method = "bilinear")
    }
    # Normalize if not already set
    if(raster::cellStats(bias, 'max') > 1 || raster::cellStats(bias, 'min') < 0 ){
      bias <- predictor_transform(bias, option = "norm")
    }
  } else { bias <- NULL }

  # Rasterize the presence estimates
  bg1 <- raster::rasterize(df[,c('x','y')] %>% sf::st_drop_geometry(),
                           background, fun = 'count', background = 0)
  bg1 <- raster::mask(bg1, background)

  assertthat::assert_that(
    is.finite(raster::cellStats(bg1,'max',na.rm = T)[1])
  )

  # Generate pseudo absence data
  if(method == "random"){
    # Now sample from all cells not occupied
    if(!is.null(bias)){
      # Get probability values for cells where no sampling has been conducted
      prob_bias <- bias[which(bg1[]==0)]
      if(any(is.na(prob_bias))) prob_bias[is.na(prob_bias)] <- 0
      abs <- sample(which(bg1[]==0), size = nrpoints, replace = TRUE, prob = prob_bias)
    } else {
      # raster::sampleStratified(bg1, nrpoints)[,1]
      abs <- sample(which(bg1[]==0), size = nrpoints, replace = TRUE)
    }
  } else if(method == "buffer"){
    assertthat::assert_that(is.numeric(buffer_distance),msg = "Buffer distance parameter not numeric!")
    bg2 <- bg1; bg2[bg2 == 0] <- NA
    # Calculate distance to all cells that are NA
    dis <- raster::distance(bg2)
    # Set values lower than XX to NA
    dis[dis <= buffer_distance] <- NA
    dis <- raster::mask(dis, background)
    dis[dis > 0] <- 0
    # Now sample from all cells not occupied
    if(!is.null(bias)){
      # Get probability values for cells where no sampling has been conducted
      prob_bias <- bias[which(dis[]==0)]
      if(any(is.na(prob_bias))) prob_bias[is.na(prob_bias)] <- 0
      abs <- sample(which(dis[]==0), size = nrpoints, replace = TRUE, prob = prob_bias)
    } else {
      # raster::sampleStratified(bg1, nrpoints)[,1]
      abs <- sample(which(dis[]==0), size = nrpoints, replace = TRUE)
    }
    rm(dis)
  } else if(method == "mcp"){
    # Idea is to draw a MCP polygon around all points and sample background points only outside of it.
    pol <- sf::st_as_sf( sf::st_convex_hull(sf::st_union(df)) )
    # Now mask out this area from the background
    inside <- settings$get("mcp_inside"); assertthat::assert_that(is.logical(inside), "MCP inside / outside parameter has to be set.")
    bg2 <- raster::mask(bg1, mask = pol, inverse = !inside)
    if(!is.null(bias)){
      # Get probability values for cells where no sampling has been conducted
      prob_bias <- bias[which(bg2[]==0)]
      if(any(is.na(prob_bias))) prob_bias[is.na(prob_bias)] <- 0
      abs <- sample(which(bg2[]==0), size = nrpoints, replace = TRUE, prob = prob_bias)
    } else {
      # raster::sampleStratified(bg1, nrpoints)[,1]
      abs <- sample(which(bg2[]==0), size = nrpoints, replace = TRUE)
    }
    rm(bg2)
  }

  # Append to the presence information.
  abs <- sf::st_as_sf(data.frame(raster::xyFromCell(bg1, abs)),
                      coords = c("x", "y"), crs = sf::st_crs(df) )
  abs$x <- sf::st_coordinates(abs)[,1]; abs$y <- sf::st_coordinates(abs)[,2]
  abs$type <- "Pseudo-absence"; abs[[field_occurrence]] <- 0
  sf::st_geometry(abs) <- attr(df, "sf_column") # Rename geom column to be the same as for df
  assertthat::assert_that( nrow(abs) > 0,
                           all(names(abs) %in% names(df)))
  # Unique to remove any duplicate values (otherwise double counted cells)
  # FIXME: Ignoring this as one might want to stress contrast to biases cells
  # abs <- unique(abs)
  # Combine with presence information and return
  out <- rbind(df, abs)
  return(out)
}

#' Aggregate count observations to a grid
#'
#' @description
#' This function aggregates provided point data to a reference grid, by,
#' depending on the type, either counting the number of observations per grid cell
#' or aggregating them via a sum.
#' @param df A [`sf`], [`data.frame`] or [`tibble`] object containing point data.
#' @param template A [`RasterLayer`] object that is aligned with the predictors.
#' @param field_occurrence A [`character`] name of the column containing the presence information (Default: \code{observed}).
#' @returns A [`sf`] object with the newly aggregated points.
#' @keywords internal
#' @noRd
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
    pres <- raster::rasterize(df, field = field_occurrence,
                              template, fun = 'sum', background = 0)

  } else {
    # Simply count them
    if(inherits(df, 'sf')) df <- df %>% sf::st_drop_geometry()
    pres <- raster::rasterize(df[,c("x","y")],
                              template, fun = 'count', background = 0)
  }
  assertthat::assert_that(
    is.Raster(pres), is.finite(raster::cellStats(pres, "max"))
  )
  if(inherits(df, 'sf')) df <- df %>% sf::st_drop_geometry()
  # Get cell ids
  ce <- raster::cellFromXY(pres, df[,c("x","y")])
  # Get new presence data
  obs <- cbind(
    data.frame(observed = raster::values(pres)[ce],
               raster::xyFromCell(pres, ce) # Center of cell
    )
  ) %>%
    # Unique to remove any duplicate values (otherwise double counted cells)
    unique()

  # Convert to sf again
  obs <- sf::st_as_sf(obs, coords = c("x", "y"), crs = sf::st_crs(df))
  obs$x <- sf::st_coordinates(obs[attr(obs, "sf_column")])[,1]
  obs$y <- sf::st_coordinates(obs[attr(obs, "sf_column")])[,2]

  # Set CRS again
  suppressWarnings(
    obs <- sf::st_set_crs(obs, value = sf::st_crs(df))
  )
  return(obs)
}

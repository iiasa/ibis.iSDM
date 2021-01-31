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

#' Logistic transformation function
#' @param x A [`numeric`] value
#' @noRd
logistic <- function(a){
  assertthat::assert_that(is.numeric(a),msg = 'Provided value not numeric')
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
#' @returns A [`list`] object with [`formula`] objects
#' @examples \dontrun{
#' formula_combinations(form)
#' }
#' @noRd

formula_combinations <- function(form, InterceptOnly = T){
  assertthat::assert_that(inherits(form,'formula'),
                          'purrr' %in% loadedNamespaces())

  # Extract variables from formula
  response = all.vars(form)[1]
  varnames = all.vars(form)[-1]

  # Remove spatial field and spde from the list
  varnames <- grep('spatial.field',varnames,value = T,invert = T)
  varnames <- grep('spde',varnames,value = T,invert = T)

  # Construct unique combinations
  varnames_comb <- 1:length(varnames) %>%
      purrr::map(~ combn(varnames, .x) %>% apply(2, list) %>% unlist(recursive = F)) %>%
    unlist(recursive = F)

  form_temp <- varnames_comb %>% purrr::map(~paste0(response, " ~ ", paste(.x, collapse = " + ")) %>% as.formula)
  if(InterceptOnly) form_temp <- form_temp %>% append(paste0(response, " ~ 1") %>% as.formula %>% list, .)

  return(form_temp)
}

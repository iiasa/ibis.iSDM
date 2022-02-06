#' @include utils.R
NULL

#' Pipe operator
#'
#' This package uses the pipe operator (`\%>\%`) to express nested code
#' as a series of imperative procedures.
#'
#' @param lhs, rhs An object and a function.
#' @seealso [magrittr::%>%()], [tee()].
#' @return An object.
#' @keywords misc
#' @examples
#' # set seed for reproducibility
#' set.seed(500)
#'
#' # generate 100 random numbers and calculate the mean
#' mean(runif(100))
#'
#' # reset the seed
#' set.seed(500)
#'
#' # repeat the previous procedure but use the pipe operator instead of nesting
#' # function calls inside each other.
#' runif(100) %>% mean()
#'
#' @name %>%
#' @rdname pipe
#' @aliases pipe
#' @importFrom magrittr %>%
#' @export
NULL

#' Central colour repository
#' @description This command contains all the colours
#' specified for use in \pkg{ibis.iSDM}.
#' @name ibis_colours
#' @examples
#' ibis_colours[['viridis_plasma']]
#' @keywords internal
#' @noRd
ibis_colours <- list(
  sdm_colour = colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'))(100),
  prob_colour = colorRampPalette(c("limegreen","springgreen4","cornflowerblue","dodgerblue4","yellow","orange","mediumvioletred","red"))(100),
  divg_bluegreen = c("#2C194C","#284577","#4B76A0","#8CA7C3","#D0DCE6","#D4E6D6","#98C39B","#5C9F61","#3E7229","#434C01"),
  divg_bluered = c("#4E193D","#44234E","#3B3264","#34487B","#376091","#4B7BA5","#6996B6","#8DADC3","#B1BEC7","#CCC1BE","#D8B7A7",
  "#D8A589","#CE8C6A","#BF724C","#A95432","#8E3821","#77231D","#661723","#5A152D","#50193B"),
  viridis_orig = c("#440154FF","#482878FF","#3E4A89FF","#31688EFF","#26828EFF","#1F9E89FF","#35B779FF","#6DCD59FF","#B4DE2CFF","#FDE725FF"),
  viridis_cividis = c("#00204DFF","#00336FFF","#39486BFF","#575C6DFF","#707173FF","#8A8779FF","#A69D75FF","#C4B56CFF","#E4CF5BFF","#FFEA46FF"),
  viridis_plasma = c("#0D0887FF","#47039FFF","#7301A8FF","#9C179EFF","#BD3786FF","#D8576BFF","#ED7953FF","#FA9E3BFF","#FDC926FF","#F0F921FF"),
  distinct_random =  c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")
)

#' Print ibis options
#'
#' @description There are a number of hidden options that can be specified for ibis.iSDM.
#' Currently supported are:
#'  * \code{'ibis.runparallel'} : [`logical`] value on whether processing should be run in parallel
#'  * \code{'ibis.nthread'} : [`numeric`] value on how many cores should be used by default
#'  * \code{'ibis.setupmessages'} : [`logical`] value indicating whether message during object creation should be shown
#'  * \code{'ibis.engines'} : Returns a [`vector`] with all valid engines.
#'  * \code{'ibis.use_future'} : [`logical`] on whether the \pkg{future} package should be used for parallel computing.
#' @return The output of \code{getOptions} for all ibis related variables
#' @keywords misc
#' @export
ibis_options <- function(){
  what <- grep('ibis',names(options()),value = TRUE)
  items <- options()[what]
  print(items)
}

#' Options to set up ibis for parallel processing with future
#'
#' @param cores A [`numeric`] number stating the number of cores to use.
#' @param strategy A [`character`] denoting the strategy to be used for future. See help of [`future`] for options.
#' (Default: \code{"multisession"}).
#' @seealso [future]
#' @keywords misc
#' @export
ibis_future <- function(cores = getOption("ibis.nthread"), strategy = "multisession") {
  assertthat::assert_that(
    is.numeric(cores),
    is.character(strategy)
  )
  strategy <- match.arg(strategy, c("sequential", "multisession", "multicore", "cluster", "remote"),
                        several.ok = FALSE)
  check_package("future")

  if(isTRUE(Sys.info()[["sysname"]] == "Windows")){
    if(strategy == "multicore") stop("Multicore is not supported ")
  }

  # Define plan
  if(strategy == "remote"){
    #TODO: See if a testing environment could be found.
    stop("TBD. Requires specific setup.")
  }
  ev <- switch(strategy,
               "sequential" = future::sequential(),
               "multisession" = future::multisession(workers = cores),
               "multicore" = future::multicore(workers = cores),
               "cluster" = future::cluster(workers = cores)
  )
  # Set up plan
  future::plan(ev)
  invisible()
}


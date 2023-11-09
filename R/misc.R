#' @include utils.R
NULL

#' Central colour repository
#' @description This command contains all the colours
#' specified for use in \pkg{ibis.iSDM}.
#' @name ibis_colours
#' @examples
#' ibis_colours[['viridis_plasma']]
#' @aliases ibis_colours
#' @keywords internal
#' @noRd
ibis_colours <- list(
  sdm_colour = colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'))(100),
  prob_colour = colorRampPalette(c("grey90","springgreen4","cornflowerblue","dodgerblue4","yellow","orange","mediumvioletred","red"))(100),
  ohsu_palette = colorRampPalette(c("white","#fbcc3f", "#56ab6c", "#5e9dcc", "#575d5f"))(100),
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
#' @description There are a number of hidden options that can be specified for
#'   ibis.iSDM. Currently supported are:
#'  * \code{'ibis.runparallel'} : [`logical`] value on whether processing should be run in parallel.
#'  * \code{'ibis.nthread'} : [`numeric`] value on how many cores should be used by default.
#'  * \code{'ibis.setupmessages'} : [`logical`] value indicating whether message during object creation should be shown.
#'  * \code{'ibis.engines'} : Returns a [`vector`] with all valid engines.
#'  * \code{'ibis.use_future'} : [`logical`] on whether the \pkg{future} package should be used for parallel computing.
#' @return The output of \code{getOptions} for all ibis related variables.
#' @keywords misc
#' @aliases ibis_options
#' @examples
#'  ibis_options()
#' @export
ibis_options <- function(){
  what <- grep('ibis',names(options()),value = TRUE)
  items <- options()[what]
  items
}

#' Install ibis dependencies
#'
#' @description Some of the dependencies (R-Packages) that ibis.iSDM relies on
#' are by intention not added to the Description of the file to keep the number
#' of mandatory dependencies small and enable the package to run even on systems
#' that might not have all libraries pre-installed.
#'
#' This function provides a convenience wrapper to install those missing
#' dependencies as needed. It furthermore checks which packages require updating
#' and updates them as needed.
#' @note INLA is handled in a special way as it is not available via cran.
#' @param deps A [`vector`] with the names of the packages to be installed
#'   (Default: \code{"ibis.dependencies"} in [`ibis_options`]).
#' @param update A [`logical`] flag of whether all (installed) packages should
#'   also be checked for updates (Default: \code{TRUE}).
#' @returns Nothing. Packages will be installed.
#' @aliases ibis_dependencies
#' @examples
#' \dontrun{
#'   # Install and update all dependencies
#'   ibis_dependencies()
#' }
#' @keywords misc
#' @export
ibis_dependencies <- function(deps = getOption("ibis.dependencies"), update = TRUE){
  assertthat::assert_that(
    is.vector(deps),
    length(deps) >= 1,
    is.logical(update)
  )
  # First check which packages are not installed and then do so.
  new.packages <- deps[!(deps %in% utils::installed.packages()[, "Package"])]
  if(length(new.packages)>0){
    if("INLA" %in% new.packages){
      if (!requireNamespace("BiocManager", quietly = TRUE))
        utils::install.packages("BiocManager")
      BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE)
      # Then install INLA
       utils::install.packages("INLA",
                         repos=c(getOption("repos"),
                                 INLA="https://inla.r-inla-download.org/R/stable"),
                         dep=TRUE)
    }
    suppressMessages(
      utils::install.packages(new.packages, dependencies = TRUE, quiet = TRUE,
                              repos = "https://cloud.r-project.org")
    )
  }

  # Addition of rangeShifter
  check_package("devtools")
  try({
    devtools::install_github("https://github.com/RangeShifter/RangeShiftR-package", ref = "main")
  })
  # Try and Add steps
  try({
    remotes::install_github("steps-dev/steps", build_vignettes = TRUE)
  })

  # Update packages if set
  if(update){
    if("INLA" %in% deps){
      # For windows
      if(length(grep("Windows", utils::osVersion, ignore.case = TRUE)) && !("INLA" %in% utils::installed.packages()[, "Package"])){
        # On windows we remove INLA and reinstall
        utils::install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
      } else if(requireNamespace("INLA", quietly = TRUE)) {
        suppressPackageStartupMessages(
          INLA::inla.upgrade(ask = FALSE)
        )
      }
    }
    # Update all the package excluding INLA
    suppressMessages(
      utils::update.packages(deps, ask = FALSE, repos = "https://cloud.r-project.org")
    )
  }
  invisible()
}

#' Options to set up ibis for parallel processing with future
#'
#' @param cores A [`numeric`] number stating the number of cores to use.
#' @param strategy A [`character`] denoting the strategy to be used for future.
#'   See help of [`future`] for options. (Default: \code{"multisession"}).
#' @seealso [future]
#' @return None
#' @aliases ibis_future
#' @examples
#' \dontrun{
#' # Starts future job
#' ibis_future(cores = 4)
#' }
#' @keywords misc
#' @export
ibis_future <- function(cores = getOption("ibis.nthread"), strategy = getOption("ibis.futurestrategy")) {
  assertthat::assert_that(
    is.numeric(cores),
    is.character(strategy)
  )
  check_package("future")
  # Check that number of cores don't exceed what is possible
  assertthat::assert_that(cores <= future::availableCores())

  strategy <- match.arg(strategy, c("sequential", "multisession", "multicore", "cluster", "remote"),
                        several.ok = FALSE)

  if(isTRUE(Sys.info()[["sysname"]] == "Windows")){
    if(strategy == "multicore") stop("Multicore is not supported on windows!")
  }

  # Define plan based on formulated strategy
  if(strategy == "remote"){
    #TODO: See if a testing environment could be found.
    stop("TBD. Requires specific setup.")
    #e.g. cl <- makeCluster(4, type = "MPI")
  } else if(strategy == "sequential") {
    future::plan(strategy = future::sequential())
  } else if(strategy == "multisession"){
    future::plan(strategy = future::multisession(workers = cores) )
  } else if(strategy == "multicore"){
    future::plan(strategy = future::multicore(workers = cores) )
  } else if(strategy == "cluster"){
    future::plan(strategy = future::cluster(workers = cores) )
  }
  # Register the doFuture adapate
  doFuture::registerDoFuture()
  invisible()
}

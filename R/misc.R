#' Central colour repository
#' @description This command contains all the colours specified for use in \pkg{ibis.iSDM}.
#'
#' @examples
#' ibis_colours[['viridis_plasma']]
#'
#' @noRd
#'
#' @keywords internal
ibis_colours <- list(
  sdm_colour = colorRampPalette(c('grey85','steelblue4','steelblue1','gold','red1','red4'))(100),
  prob_colour = colorRampPalette(c("grey85","springgreen4","cornflowerblue","dodgerblue4","yellow","orange","mediumvioletred","red"))(100),
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
#' @description There are a number of hidden options that can be specified for ibis.iSDM.
#' Currently supported are:
#' * \code{'ibis.runparallel'} : [`logical`] value on whether processing should be run in parallel.
#' * \code{'ibis.nthread'} : [`numeric`] value on how many cores should be used by default.
#' * \code{'ibis.setupmessages'} : [`logical`] value indicating whether message during object creation should be shown (Default: \code{NULL}).
#' * \code{'ibis.engines'} : Returns a [`vector`] with all valid engines.
#' * \code{'ibis.use_future'} : [`logical`] on whether the \pkg{future} package should be used for parallel computing.
#'
#' @return The output of \code{getOptions} for all ibis related variables.
#'
#' @keywords misc
#'
#' @examples
#' ibis_options()
#'
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
#'
#' @param deps A [`vector`] with the names of the packages to be installed
#' (Default: \code{"ibis.dependencies"} in [`ibis_options`]).
#' @param update A [`logical`] flag of whether all (installed) packages should
#' also be checked for updates (Default: \code{TRUE}).
#'
#' @note INLA is handled in a special way as it is not available via cran.
#'
#' @returns Nothing. Packages will be installed.
#'
#' @keywords misc
#'
#' @examples
#' \dontrun{
#'   # Install and update all dependencies
#'   ibis_dependencies()
#' }
#'
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

#' Set the parallel processing flag to TRUE
#' @description
#' Small helper function to enable parallel processing. If set
#' to \code{TRUE}, then parallel inference (if supported by engines) and projection is
#' enabled across the package.
#' For enabling prediction support beyond sequential prediction see the [`ibis_future`] function.
#'
#' @return Invisible
#' @seealso [future], [ibis_future]
#' @keywords misc
#' @export
ibis_enable_parallel <- function(){
  options('ibis.runparallel' = TRUE)
  if(getOption('ibis.nthread')<2){
    myLog('[Setup]','yellow','Parallelization enabled but less than 2 nodes specified!')
  }
  invisible()
}

#' Set the number of threads for parallel processing.
#' @description
#' Small helper function to respecify the strategy for parallel processing (Default: \code{'sequential'}).
#' @details
#' Currently supported strategies are:
#'
#' * \code{"sequential"} = Resolves futures sequentially in the current R process (Package default).
#' * \code{"multisession"} = Resolves futures asynchronously across \code{'cores'} sessions.
#' * \code{"multicore"} = Resolves futures asynchronously across on forked processes. Only works on UNIX systems!
#' * \code{"cluster"} = Resolves futures asynchronously in sessions on this or more machines.
#' * \code{"slurm"} = To be implemented: Slurm linkage via batchtools.
#' @param strategy A [`character`] with the strategy.
#' @return Invisible
#' @seealso [future], [ibis_future]
#' @keywords misc
#' @export
ibis_set_strategy <- function(strategy = "sequential"){
  assertthat::assert_that(is.character(strategy))

  strategy <- match.arg(strategy, c("sequential", "multisession", "multicore", "cluster", "slurm"),
                        several.ok = FALSE)
  options('ibis.futurestrategy' = strategy)
  invisible()
}

#' Set the threads for parallel processing.
#' @description
#' Small helper function to respecify the number of threads for parallel processing.
#' @param threads A [`numeric`] greater thna \code{0}.
#' @return Invisible
#' @seealso [future], [ibis_future_run]
#' @keywords misc
#' @export
ibis_set_threads <- function(threads = 2){
  assertthat::assert_that(is.numeric(threads),
                          threads >0)
  options('ibis.nthread' = threads)
  invisible()
}

#' Internal function to enable  (a)synchronous parallel processing
#'
#' @description
#' This function checks if parallel processing can be set up and enables it.
#' **Ideally this is done by the user for more control!**
#' In the package parallelization is usually only used for predictions and projections,
#' but not for inference in which case parallel inference should be handled by the engine.
#' @details
#' Currently supported strategies are:
#'
#' * \code{"sequential"} = Resolves futures sequentially in the current R process (Package default).
#' * \code{"multisession"} = Resolves futures asynchronously across \code{'cores'} sessions.
#' * \code{"multicore"} = Resolves futures asynchronously across on forked processes. Only works on UNIX systems!
#' * \code{"cluster"} = Resolves futures asynchronously in sessions on this or more machines.
#' * \code{"slurm"} = To be implemented: Slurm linkage via batchtools.
#' @note
#' The \code{'plan'} set by [future] exists after the function has been executed.
#'
#' If the aim is to parallize across many species, this is better done in a scripted solution.
#' Make sure not to parallize predictions within existing clusters to avoid out-of-memory
#' issues.
#'
#' @param plan_exists A [`logical`] check on whether an existing [`future`] plan exists (Default: \code{FALSE}).
#' @param cores A [`numeric`] number stating the number of cores to use.
#' @param strategy A [`character`] denoting the strategy to be used for future.
#' See help of [`future`] for options. (Default: \code{"multisession"}).
#' @param workers An optional list of remote machines or workers, e.g. \code{"c(remote.server.org)"}.
#' Alternatively a \code{"cluster"} object can be provided.
#'
#' @return Invisible
#'
#' @seealso [future]
#' @keywords misc
#'
#' @examples
#' \dontrun{
#' # Starts future job. F in this case is a prediction function.
#' ibis_future(cores = 4, strategy = "multisession")
#' }
#'
#' @export
ibis_future <- function(plan_exists = FALSE,
                        cores = getOption("ibis.nthread",default = 2),
                        strategy = getOption("ibis.futurestrategy"),
                        workers = NULL
                        ) {
  assertthat::assert_that(
    is.logical(plan_exists),
    is.numeric(cores),
    is.character(strategy),
    is.null(workers) || is.vector(workers) || (inherits(workers, "ClusterFuture"))
  )

  # Check if plan exists, if not specify new plan
  if(!plan_exists){
    # Check that number of cores don't exceed what is possible
    assertthat::assert_that(cores <= parallelly::availableCores()[[1]])

    strategy <- match.arg(strategy, c("sequential", "multisession", "multicore", "cluster", "slurm"),
                          several.ok = FALSE)

    # Check if parallel processing is enabled
    assertthat::assert_that(
      getOption("ibis.runparallel", default = FALSE),
      msg = "Parallel processing not enabled. Run 'ibis_enable_parallel()' !"
    )

    # Check that more 1 connection is available
    if(strategy != "sequential"){
      assertthat::assert_that(parallelly::availableConnections()>1,
                              msg = "No further connections are available for parallel processing!")
    }

    # isTRUE(Sys.info()[["sysname"]] == "Windows")
    if(strategy == "multicore" && !parallelly::supportsMulticore()){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Setup]','yellow','Parallization multicore not supported om windows. Changing to multisession.')
      strategy <- "multisession"
    }

    if(is.null(cores)) cores <- 4 # Arbitrary nr of cores if somehow not found
    if(is.null(workers)) workers <- cores # If no workers are found, use the cores

    # --- #
    # Define plan based on formulated strategy
    if(strategy == "slurm"){
      #TODO: See if a testing environment could be found.
      cli::cli_abort("Not yet implemented")
      #e.g. cl <- makeCluster(4, type = "MPI")
    } else if(strategy == "sequential") {
      future::plan(strategy = "sequential")
    } else if(strategy == "multisession"){
      future::plan(strategy = "multisession", workers = cores)
    } else if(strategy == "multicore"){
      future::plan(strategy = "multicore", workers = cores)
    } else if(strategy == "cluster"){
      future::plan(strategy = "cluster", workers = cores)
    }
  }
  if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Setup]','green','Specified parallel processing plan with strategy: ', strategy)
  invisible()
}

#' Internal helper function to split a data.frame into chucnks
#'
#' @param X A [`data.frame`] or [`matrix`] object to be split.
#' @param N A [`numeric`] with the number of chunks smaller than \code{`nrow(X)`}.
#' @param cores A [`numeric`] with the number of processing cores. Used when \code{N} is \code{NULL}.
#' @param index_only [`logical`] on whether only the indices or split X as list is returnsed (Default: \code{FALSE}).
#' @keywords internal
#' @returns A [`list`] object.
#' @examples
#' # Chunck example data into 4 subsets
#' chunk_data(datasets::airquality, N = 4)
#' @noRd
chunk_data <- function(X, N = NULL, cores = parallel::detectCores(), index_only = FALSE){
  assertthat::assert_that(
    is.data.frame(X) || is.matrix(X),
    nrow(X) > 1,
    is.numeric(cores),
    is.null(N) || is.numeric(N),
    is.logical(index_only)
  )

  n_vars <- nrow(X)
  # Use cores as N otherwise
  if(is.null(N)) N <- cores
  chunk_size <- ceiling(n_vars / N)
  n_chunks <- ceiling(n_vars / chunk_size)

  if(index_only){
    chunk_list <- list()
  } else {
    chunk_list <- vector(length = n_chunks, mode = "list")
  }
  for (i in seq_len(n_chunks)) {
    if ((chunk_size * (i - 1) + 1) <= n_vars) {
      chunk <- (chunk_size * (i - 1) + 1):(min(c(chunk_size *
                                                   i, n_vars)))
      if(index_only){
        o <- chunk
      } else {
        o <- X[chunk, ] |> as.data.frame()
        colnames(o) <- colnames(X)
      }
      chunk_list[[i]] <- o
      rm(o)
    }
  }
  assertthat::assert_that(is.list(chunk_list))
  if(!index_only){
    assertthat::assert_that(sum(sapply(chunk_list, nrow)) == nrow(X),
                            msg = "Something went wrong with the data chunking...")
  }
  return(chunk_list)
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
#' @param export_packages A [`vector`] with packages to export for use on parallel
#' nodes (Default: \code{NULL}).
#' @param ... Any other parameter passed on.
#'
#' @details By default, the [parallel] package is used for parallel computation,
#' however an option exists to use the [future] package instead.
#'
#' @keywords utils
#'
#' @examples
#' \dontrun{
#'  run_parallel(list, mean, cores = 4)
#' }
#'
#' @export
run_parallel <- function(X, FUN, cores = 1, approach = "future", export_packages = NULL, ...) {
  assertthat::assert_that(
    is.list(X) || is.data.frame(X) || is.matrix(X),
    is.function(FUN),
    is.numeric(cores),
    is.null(export_packages) || is.character(export_packages)
  )
  message("The run_parallel function is likely deprecated and is only kept for reference...")

  # Match approach
  approach <- match.arg(approach, c("parallel", "future"), several.ok = FALSE)

  # Collect dots
  dots <- list(...)

  if(!is.list(X)){
    # Convert input object to a list of split parameters
    X <- chunk_data(X, cores = cores)
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
      check_package("future")
      ibis_future()
    }
  }
  # If input data was not a list, combine again
  if(input_type != "list" && is.list(out)){
    out <- do.call(rbind, out)
  }
  return( out )
}

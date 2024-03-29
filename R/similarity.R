#' @include class-biodiversitydistribution.R class-distributionmodel.R
NULL

#' Calculate environmental similarity of reference datasets to predictors.
#'
#' @description Calculate the environmental similarity of the provided
#' covariates with respect to a reference dataset. Currently supported is
#' Multivariate Environmental Similarity index and the multivariate combination
#' novelty index (NT2) based on the Mahalanobis divergence (see references).
#'
#' @param obj A [`BiodiversityDistribution`], [`DistributionModel`] or alternatively
#' a [`SpatRaster`] object.
#' @param ref A [`BiodiversityDistribution`], [`DistributionModel`] or alternatively
#' a [`data.frame`] with extracted values (corresponding to those given in `obj`).
#' @param ref_type A [`character`] specifying the type of biodiversity to use
#' when obj is a [`BiodiversityDistribution`].
#' @param method A specifc method for similarity calculation. Currently supported:
#' \code{'mess'}, \code{'nt'}.
#' @param predictor_names An optional [`character`] specifying the covariates to
#' be used (Default: \code{NULL}).
#' @param full should similarity values be returned for all variables (Default:\code{FALSE})?
#' @param plot Should the result be plotted? Otherwise return the output list (Default: \code{TRUE}).
#' @param ... other options (Non specified).
#'
#' @details [`similarity`] implements the MESS algorithm described in Appendix S3
#' of Elith et al. (2010) as well as the Mahalanobis dissimilarity described in
#' Mesgaran et al. (2014).
#'
#' @returns This function returns a list containing:
#' * `similarity`: A `SpatRaster` object with multiple layers giving the environmental
#' similarities for each variable in `x` (only included when \code{"full=TRUE"});
#'  * `mis`: a `SpatRaster` layer giving the minimum similarity value
#' across all variables for each location (i.e. the MESS);
#'  * `exip`: a `SpatRaster` layer indicating whether any model would interpolate
#' or extrapolate to this location based on environmental surface;
#'  * `mod`: a factor `SpatRaster` layer indicating which variable was most
#' dissimilar to its reference range (i.e. the MoD map, Elith et al. 2010); and
#'  * `mos`: a factor `SpatRaster` layer indicating which variable was most
#' similar to its reference range.
#'
#' @references
#' * Elith, J., Kearney, M., and Phillips, S. (2010) "The art of modelling range-shifting
#' species". _Methods in Ecology and Evolution_, 1: 330-342. https://doi.org/10.1111/j.2041-210X.2010.00036.x
#' * Mesgaran, M.B., Cousens, R.D. and Webber, B.L. (2014) "Here be dragons: a tool
#' for quantifying novelty due to covariate range and correlation change when
#' projecting species distribution models". _Diversity and Distributions_, 20: 1147-1159.
#' https://doi.org/10.1111/ddi.12209
#'
#' @seealso dismo R-package.
#' @keywords mess mahalanobis similarity environment
#'
#' @examples
#'  \dontrun{
#' plot(
#'   similarity(x) # Where x is a distribution or Raster object
#' )
#'  }
#'
#' @name similarity
NULL

#' @rdname similarity
#' @export
methods::setGeneric(
  "similarity",
  signature = methods::signature("obj"),
  function(obj, ref, ref_type = 'poipo',
           method = 'mess', predictor_names = NULL, full = FALSE, plot = TRUE, ...) standardGeneric("similarity"))

#' Similarity of used predictors from a trained distribution model
#' @rdname similarity
methods::setMethod(
  "similarity",
  methods::signature(obj = "BiodiversityDistribution"),
  function(obj, ref_type = 'poipo', method = 'mess', predictor_names = NULL, full = FALSE, plot = TRUE, ...) {
    assertthat::assert_that(inherits(obj, "BiodiversityDistribution"),
                            is.character(ref_type),
                            ref_type %in% c('poipo','poipa','polpo','polpa'),
                            is.character(method),
                            is.null(predictor_names) || is.character(predictor_names)
    )
    # ref_type = 'poipo'; method = 'mess'; predictor_names = NULL; full = FALSE
    # Check that data and predictors are there
    assertthat::assert_that(!is.Waiver(obj$biodiversity),!is.Waiver(obj$predictors))
    assertthat::assert_that(is.null(predictor_names) || all(predictor_names %in% obj$get_predictor_names()) )
    # Match to correct spelling mistakes
    method <- match.arg(tolower(method), c('mess','nt'), several.ok = FALSE)

    # Get biodiversity data
    bid <- obj$biodiversity$get_id_byType(ref_type)
    if(length(bid)==0){
      # Try and find an alternative layer
      if(obj$biodiversity$length() == 1){
        bid <- obj$biodiversity$get_types()
      } else stop("No biodiversity layer found of the given type. Please check biodversity data!")
    }

    # Get covariates
    covs <- obj$predictors$get_data()
    assertthat::assert_that(is.Raster(covs), msg = "No Raster layer found!")
    # Extract covariates for reference data
    ref <- get_rastervalue(coords = obj$biodiversity$get_coordinates(names(bid)),
                           env = covs,
                           rm.na = FALSE)
    ref <- ref[,names(covs)]

    # Subset if necessary
    if(!is.null(predictor_names)){ covs <- covs[[predictor_names]]; ref <- ref[,predictor_names]}

    assertthat::assert_that( terra::nlyr(covs) == ncol(ref))

    # Run mess function
    if(method == 'mess'){
      out <- .mess(covs = covs,
                   ref  = ref,
                   full = full)
      # Calculate interpolation/extrapolated
      rip <- terra::classify(out$mis,
                             c( terra::global(out$mis,'min', na.rm = TRUE)[,1], 0,
                                   terra::global(out$mis,'max', na.rm = TRUE)[,1]))
      rip <- terra::as.factor(rip)
      for(i in 1:terra::nlyr(rip)){
        ca <- data.frame(ID = levels(rip[[i]])[[1]][,1])
        ca[names(rip[[i]])] <- c('Extrapolation','Interpolation')
        levels(rip[[i]]) <- ca
      }
      out$exip <- rip;rm(rip)

      # Relabel most important
      out$mod <- terra::as.factor(out$mod)
      levels(out$mod) <- data.frame(ID = levels(out$mod)[[1]][,1],
                                    variable = names(covs))
      # Reduce to SpatRaster if not one already
      if(!is.Raster(out)){
        out <- Reduce(c, out)
        names(out) <- c("mis", "mod", "mos", "exip")
      }

    } else if(method == 'nt') {
      out <- .nt12(prodat = covs,
                   refdat = ref)
    }

    # If plot is specified, make figures. Otherwise return the list of rasters
    if(plot){
      if(method == 'mess'){
        par.ori <- graphics::par(no.readonly = TRUE)
        graphics::par(mfrow=c(2,2))
        terra::plot(exp(out$mis),col = ibis_colours[['viridis_plasma']],main = paste0('Similarity surface\n (method: ',method,', exp. transformed)'))
        terra::plot(out$exip,col = ibis_colours[['distinct_random']][1:2],main = paste0('Extrapolated vs interpolated conditions'))
        terra::plot(out$mod,col = ibis_colours[['distinct_random']][1:length(unique(out$mod)[,1])], main = paste0('Most dissimilar from reference'))
        terra::plot(out$mos,col = ibis_colours[['distinct_random']][length(ibis_colours[['distinct_random']]):(length(ibis_colours[['distinct_random']])-length(unique(out$mos)[,1]))], main = paste0('Most similar to reference'))
        graphics::par(par.ori)
      } else if(method == 'nt'){
        par.ori <- graphics::par(no.readonly = TRUE)
        graphics::par(mfrow=c(3,1))
        terra::plot(exp(out$NT1),col = ibis_colours[['viridis_plasma']],main = paste0('Univariate extrapolation\n(exp. transformed)'))
        terra::plot(log(out$NT2),col = ibis_colours[['viridis_orig']],main = paste0('Non-analogous dissimilarity\n(log. transformed)'))
        terra::plot(out$novel,col = ibis_colours[['distinct_random']][4:7],main = paste0('Novel conditions (method: ',method,')'))
        graphics::par(par.ori)
      }
    } else {
      return( out )
    }
  }
)

#' Similarity of used predictors by providing a SpatRaster directly
#' @rdname similarity
methods::setMethod(
  "similarity",
  methods::signature(obj = "SpatRaster"),
  function(obj, ref, method = 'mess', full = FALSE, plot = TRUE, ...) {
    assertthat::assert_that(!missing(ref),msg = 'Provide a sf object of reference sites')
    assertthat::assert_that(inherits(obj, "SpatRaster"),
                            terra::nlyr(obj)>=1,
                            inherits(ref, 'sf'),
                            is.character(method)
    )
    # Check that points are of same projection as raster
    assertthat::assert_that(sf::st_crs(obj) == sf::st_crs(ref))

    # Match to correct spelling mistakes
    method <- match.arg(tolower(method), c('mess','nt2'), several.ok = FALSE)

    # Extract values for each provided value
    ex <- terra::extract(x = obj,
                          y = ref,
                          df = TRUE)
    # Subset to variables in obj and remove missing rows
    ex <- subset.data.frame(ex, select = names(obj))
    ex <- subset.data.frame(ex, stats::complete.cases(ex))

    if(method == 'mess'){
      out <- .mess(covs = obj,
                   ref  = ex,
                   full = full)
      # Calculate interpolation/extrapolated
      rip <- terra::classify(out$mis,
                           c(terra::global(out$mis,'min', na.rm = TRUE)[,1], 0,
                             terra::global(out$mis,'max', na.rm = TRUE)[,1]),
                           na.rm = TRUE)
      rip <- terra::as.factor(rip)
      for(i in 1:terra::nlyr(rip)){
        ca <- data.frame(ID = levels(rip[[i]])[[1]][,1])
        ca[names(rip[[i]])] <- c('Extrapolation','Interpolation')
        levels(rip[[i]]) <- ca
      }
      out$exip <- rip;rm(rip)

      # Relabel most important
      out$mod <- terra::as.factor(out$mod)
      levels(out$mod) <- data.frame(ID = levels(out$mod)[[1]][,1],
                                    variable = names(obj))

    } else {
      stop('Not yet implemented!')
    }

    # If plot is specified, make figures. Otherwise return the list of rasters
    if(plot){
      if(method == 'mess'){
        par.ori <- graphics::par(no.readonly = TRUE)
        graphics::par(mfrow=c(2,2))
        terra::plot(out$mis, col = ibis_colours[['viridis_plasma']],
                    main = paste0('Similarity surface (method: ',method,')'))
        terra::plot(out$exip,col = ibis_colours[['distinct_random']][1:2],
                    main = paste0('Extrapolated vs interpolated conditions'))
        terra::plot(out$mod,
                    col = ibis_colours[['distinct_random']][1:length(unique(out$mod))],
                    main = paste0('Most dissimilar from reference'))
        terra::plot(out$mos,
                    col = ibis_colours[['distinct_random']][length(ibis_colours[['distinct_random']]):(length(ibis_colours[['distinct_random']])-length(unique(out$mos)))],
                    main = paste0('Most similar to reference'))
        graphics::par(par.ori)
      } else if(method == 'nt'){
        par.ori <- graphics::par(no.readonly = TRUE)
        graphics::par(mfrow=c(1,3))
        terra::plot(out$NT1,col = ibis_colours[['viridis_plasma']],
                    main = paste0('Univariate extrapolation'))
        terra::plot(out$NT2,col = ibis_colours[['viridis_orig']],
                    main = paste0('Non-analogous dissimilarity'))
        terra::plot(out$novel,col = ibis_colours[['distinct_random']][1:3],
                    main = paste0('Novel conditions (method: ',method,')'))
        graphics::par(par.ori)
      }
    } else {
      return( out )
    }
  }
)

#' Function to calculate the multivariate combination novelty index (NT2)
#'
#' @description NT1 ranges from infinite negative values to zero where zero
#' indicates no extrapolation beyond the univariate coverage of reference data"
#' (Mesgaran et al. 2014).
#'
#' "NT2 can range from zero up to unbounded positive values. NT2 values ranging
#' from zero to one indicate similarity (in terms of both univariate range and
#' multivariate combination), with values closer to zero being more similar.
#' Values larger than one are indicative of novel combinations" (Mesgaran et al.
#' 2014).
#'
#' @param prodat A [`SpatRaster`]. The projected values. The layer names must
#' match the column names of \code{refdat}.
#' @param refdat A numerical [`matrix`] or [`data.frame`]. The reference values
#' of variables organized in columns.
#'
#' @note: The code is adapted from Bell & Schlaepfer 2015
#' (available at \url{https://github.com/bellland/SDM.Virtual.Species_Bell.Schlaepfer})
#' which was based on a comment by Matthew Bayly made at
#' \url{https://pvanb.wordpress.com/2014/05/13/a-new-method-and-tool-exdet-to-evaluate-novelty-environmental-conditions/}.
#'
#' @references
#'  * Mesgaran, M. B., R. D. Cousens, B. L. Webber, and J. Franklin.
#' 2014. Here be dragons: a tool for quantifying novelty due to covariate range
#' and correlation change when projecting species distribution models. Diversity
#' and Distributions 20:1147-1159.
#'
#' @noRd
#'
#' @keywords internal
.nt12 <- function(prodat, refdat){
  check_package("matrixStats")
  assertthat::assert_that(
    is.data.frame(refdat) || is.matrix(refdat)
  )
  # If not matching, check if this can be corrected
  if(terra::nlyr(prodat) != ncol(refdat)){
    refdat <- subset(refdat, select = names(prodat))
  }

  # Input checks
  assertthat::assert_that(is.Raster(prodat),
                          terra::nlyr(prodat) == ncol(refdat))
  # Make a background layer for filling
  bg <- emptyraster(prodat)
  # Now convert both to matrix
  prodat <- terra::as.matrix(prodat)
  refdat <- as.matrix(refdat)
  # Further checks
  assertthat::assert_that(identical(colnames(refdat), colnames(prodat)),
                          is.data.frame(prodat) || is.matrix(prodat),
                          is.data.frame(prodat) || is.matrix(prodat))

  # First calculate univariate novelty, e.g. NT1
  # Get ranges of variables and multiply
  range_ref <- t(matrixStats::colRanges(as.matrix(refdat), na.rm = TRUE))
  diffs_ref <- matrixStats::colDiffs(range_ref)
  # Remove those with 0 range, e.g. singular values for the reference
  range_ref <- range_ref[,which(diffs_ref!=0)]
  refdat <- refdat[,which(diffs_ref!=0)]
  prodat <- prodat[,which(diffs_ref!=0)]
  diffs_ref <- diffs_ref[,which(diffs_ref!=0)]

  range_ref_arr <- array(range_ref, dim = c(dim(range_ref), nrow(prodat)),
                         dimnames = list(c("min", "max"), colnames(refdat), NULL))

  diffs_ref_arr <- matrix(diffs_ref, nrow = nrow(prodat), ncol = ncol(prodat),
                          byrow = TRUE)

  # Make empty raster with 3 dimensions
  iud <- array(0, dim = c(dim(prodat), 3))
  iud[ , , 2] <- prodat - t(range_ref_arr["min", ,])
  iud[ , , 3] <- t(range_ref_arr["max", ,]) - prodat

  # Univariate novelty
  # NT1 ranges from infinite negative values to zero where zero indicates no extrapolation
  # beyond the univariate coverage of reference data.
  UDs <- apply(iud, 1:2, min) / diffs_ref_arr
  nt1 <- emptyraster(bg)
  nt1[] <- rowSums(UDs)

  # --- #
  # Multivariate combination novelty index (NT2)
  # Calculate the center of reference data: average and covariance matrix
  ref_av  <- colMeans(refdat, na.rm = TRUE)
  ref_cov <- stats::var(refdat, na.rm = TRUE)

  # Mahalanobis distance of reference data to center of reference data
  mah_ref <- stats::mahalanobis(x = refdat, center = ref_av, cov = ref_cov,tol=1e-20)
  # Mahalanobis distance of projected data to center of reference data
  mah_pro <- stats::mahalanobis(x = prodat, center = ref_av, cov = ref_cov,tol=1e-20)
  # Correction when mah_pro is negative (numerical instability. Correct to 0)
  if(min(mah_pro,na.rm = T) < 0) mah_pro[which(!is.na(mah_pro) & mah_pro < 0)] <- 0

  # Ratio
  mah_max <- max(mah_ref[is.finite(mah_ref)])
  nt2 <- emptyraster(bg)
  nt2[] <- (mah_pro / mah_max)

  # Calculate most dissimilar value (MOD)
  # FIXME: Implement when have time https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.12209
  # mod <- emptyraster(bg)
  # For any point, the MIC is the covariate which produces the highest ICp value.
  # icp <- 100 * (mah_ref - icp)
  # matrixStats::rowMaxs(icp)

  # --- #
  # Calculate areas outside the univariate range of combinations and
  # non-analogous novel combinations
  nt_novel <- terra::init(bg, NA)
  # First areas areas in the projection space with at least one covariate
  # outside the univariate range of reference data
  if(terra::hasValues(nt1)) o_low <- nt1 < 0 else o_low <- terra::init(nt1, 0)
  # Next areas with NT2 ranging from 0 to 1 that are similar to the reference data
  o_mid <- nt2 %in% c(0,1)
  # non-analogous covariate combinations
  o_high <- nt2 > 1

  nt_novel[o_low == 0] <- 0
  nt_novel[o_low == 1] <- 1
  nt_novel[o_mid == 1] <- 2
  nt_novel[o_high == 1] <- 3
  nt_novel <- terra::as.factor(nt_novel)
  levels(nt_novel) <- data.frame(ID = c(0,1,2,3),
                            what = c('Reference','Within reference',
                                     'Outside reference','Novel combinations'))

  # Create output stack
  out <- c(nt1, nt2, nt_novel)
  names(out) <- c('NT1','NT2','novel')
  return(out)
}

#' Function to calculate Multivariate Environmental Similarity index
#'
#' @description Internal function to calculate the MESS
#'
#' @param covs A [`SpatRaster`] with the covariates.
#' @param ref A [`data.frame`] with the covariates for the reference values.
#' @param full A [`logical`] indication whether the full extent be calculated.
#'
#' @returns Á [`SpatRaster`] object.
#'
#' @noRd
#'
#' @keywords internal
.mess <- function(covs, ref, full=FALSE) {
  assertthat::assert_that(
    is.data.frame(ref) || is.matrix(ref),
    nrow(ref)>0,
    is.logical(full)
  )
  # If not matching, check if this can be corrected
  if(terra::nlyr(covs) != ncol(ref)){
    ref <- subset(ref, select = names(covs))
  }
  assertthat::assert_that(terra::nlyr(covs) == ncol(ref))

  # Convert to data.frame
  if(!is.data.frame(ref)) {
    ref <- as.data.frame(ref, na.rm = FALSE)
  }
  # Make dummy template rasters
  if(is.Raster(covs)) {
    r <- TRUE
    if(isTRUE(full)) {
      out <- terra::init(covs, NA)
    } else {
      out <- emptyraster(covs)
    }
  } else r <- FALSE
  ref <- stats::na.omit(ref)   # Remove NAs
  if(!is.data.frame(covs)) {
    covs <- as.data.frame(covs, na.rm = FALSE)
  }
  # Calculate dimensions (range)
  if(is.null(dim(ref))) {
    rng <- as.data.frame(range(ref, na.rm=TRUE))
  } else {
    rng <- as.data.frame(apply(ref, 2, range, na.rm=TRUE))
  }
  # remove variables where max-min is 0
  rng <- rng[which(apply(rng, 2, diff)>0)]
  covs <- covs[,names(rng)]
  ref <- ref[,names(rng)]

  # Find intervals within ranges
  pct_less <- mapply(function(x, ref) {
    findInterval(x, sort(ref))/length(ref)
  }, covs, ref, SIMPLIFY=FALSE)
  # Calculate similarity surface
  sim <- mapply(function(f, rng, p) {
    ifelse(f==0, (p-rng[1])/diff(rng)*100,
           ifelse(f > 0 & f <= 0.5, f*200,
                  ifelse(f > 0.5 & f < 1, (1-f)*200,
                         (rng[2]-p)/diff(rng)*100)))
  }, pct_less, rng, covs)

  min_sim <- if(is.matrix(sim)) apply(sim, 1, min) else(min(sim))

  # Get minimum similarity and most (dis)similiar values
  mins <- apply(sim, 1, which.min)
  most_dissimilar_vec <- unlist(ifelse(lengths(mins)==0, NA, mins))
  maxs <- apply(sim, 1, which.max)
  most_similar_vec <- unlist(ifelse(lengths(maxs)==0, NA, maxs))

  if(isTRUE(r)) {
    # Calculate most dissimilar surface
    most_dissimilar <- emptyraster(out)
    most_dissimilar[] <- most_dissimilar_vec
    most_dissimilar <- as.factor(most_dissimilar)
    rat <- levels(most_dissimilar)[[1]]
    rat$varname <- colnames(sim)[rat$ID]
    levels(most_dissimilar) <- rat

    # Calculate most similar surface
    most_similar <- terra::rast(out)
    most_similar[] <- most_similar_vec
    most_similar <- as.factor(most_similar)
    rat <- levels(most_similar)[[1]]
    rat$varname <- colnames(sim)[rat$ID]
    levels(most_similar) <- rat

    # Fill template rasters
    out_min <- terra::rast(out)
    out_min[] <- min_sim
    if(isTRUE(full)) {
      out[] <- sim
      list(similarity=out, mis=out_min, mod=most_dissimilar,
           mos=most_similar)
    } else list(mis=out_min, mod=most_dissimilar, mos=most_similar)
  } else {
    if(isTRUE(full)) {
      list(similarity=sim, mis=min_sim,
           mod=most_dissimilar_vec, mos=most_similar_vec)
    } else list(mis=min_sim, mod=most_dissimilar_vec,
                mos=most_similar_vec)
  }
}

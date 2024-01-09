#' @title ibis.iSDM
#'
#' @description
#' Integrated framework of modelling the distribution of species and ecosystems in
#' a suitability framing. This package allows the estimation of integrated species
#' distribution models (iSDM) based on several sources of evidence and provided
#' presence-only and presence-absence datasets. It makes heavy use of point-process
#' models for estimating habitat suitability and allows to include spatial latent
#' effects and priors in the estimation. To do so 'ibis.iSDM' supports a number
#' of engines for Bayesian and more non-parametric machine learning estimation.
#' Further, the 'ibis.iSDM' is specifically customized to support spatial-temporal
#' projections of habitat suitability into the future.
#'
#' @name ibis.iSDM
#' @docType package
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom foreach %do% %dopar%
#' @importFrom stats effects
## usethis namespace: end
NULL

globalVariables(c("background", "band", "bi_class", "bias",
                  "change", "cid", "cell", "cluster",
                  "data",
                  "form2",
                  "geometry",
                  "id", "included", "i",
                  "km", "vt",
                  "limit", "lower", "layer",
                  "median", "model", "mad",
                  "name",
                  "observed", "oversampled",
                  "partial_effect", "polpo", "predictor", "predictors",
                  "q05", "q50", "q95",
                  "ras", "region.poly", "rug",
                  "s", "state", "suitability", "sd",
                  "tothin", "type", "time",
                  "upper",
                  "var1", "var2", "value", "variable",
                  "x", "y", "z", "zone",
                  "."))

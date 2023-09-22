#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom foreach %do% %dopar%
#' @importFrom methods as
#' @importFrom terra mask
#' @importFrom stats effects
#' @importFrom stats residuals
#' @importFrom stats complete.cases
#' @importFrom stats mad
#' @importFrom stats sd
#' @importFrom stats terms.formula
#' @importFrom stats update.formula
#' @importFrom utils install.packages
#' @importFrom graphics par
## usethis namespace: end
NULL

globalVariables(c("background", "band", "bi_class", "bias",
                  "change", "cid",
                  "data",
                  "form2",
                  "geometry",
                  "id", "included",
                  "km",
                  "limit", "lower",
                  "median", "model",
                  "name",
                  "observed", "oversampled",
                  "partial_effect", "polpo", "predictor", "predictors",
                  "q05", "q50", "q95",
                  "ras", "region.poly",
                  "s", "state", "suitability",
                  "tothin", "type", "time",
                  "upper",
                  "layer", "rug",
                  "var1", "var2", "value", "variable",
                  "x", "y", "z",
                  "."))

# Validation of a fitted distribution object

This function conducts a model evaluation based on either on the fitted
point data or any supplied independent. **Currently only supporting
point datasets. For validation of integrated models more work is
needed.**

## Usage

``` r
validate(
  mod,
  method = "continuous",
  layer = "mean",
  point = NULL,
  point_column = "observed",
  field_occurrence = NULL,
  ...
)

# S4 method for class 'ANY'
validate(
  mod,
  method = "continuous",
  layer = "mean",
  point = NULL,
  point_column = "observed",
  field_occurrence = NULL,
  ...
)

# S4 method for class 'SpatRaster'
validate(
  mod,
  method = "continuous",
  layer = NULL,
  point = NULL,
  point_column = "observed",
  field_occurrence = NULL,
  ...
)
```

## Arguments

- mod:

  A fitted
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
  object with set predictors. Alternatively one can also provide
  directly a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  however in this case the `point` layer also needs to be provided.

- method:

  Should the validation be conducted on the continious prediction or a
  (previously calculated) thresholded layer in binary format? Note that
  depending on the method different metrics can be computed. See
  Details.

- layer:

  In case multiple layers exist, which one to use? (Default: `'mean'`).

- point:

  A [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) object
  with type `POINT` or `MULTIPOINT`.

- point_column:

  A [`character`](https://rdrr.io/r/base/character.html) vector with the
  name of the column containing the independent observations. (Default:
  `'observed'`).

- field_occurrence:

  (Deprecated) A [`character`](https://rdrr.io/r/base/character.html)
  field pointing to the name of the independent observations. Identical
  to `"point_column"`

- ...:

  Other parameters that are passed on. Currently unused.

## Value

Return a tidy
[`tibble::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
with validation results.

## Details

The `'validate'` function calculates different validation metrics
depending on the output type.

The output metrics for each type are defined as follows: (where TP
stands for true positive, TN for true negative, FP the false positive
and FN the false negative) **Continuous:**

- `'n'` = Number of observations.

- `'rmse'` = Root Mean Square Error (RMSE): \$\$RMSE = \sqrt{\frac{1}{N}
  \sum\_{i=1}^{N} (\hat{y}\_{i} - y\_{i})^2}\$\$

- `'mae'` = Mean Absolute Error (MAE): \$\$MAE = \frac{\sum\_{i=1}^{N}
  \|y\_{i} - x\_{i}\|}{N}\$\$

- `'logloss'` = Log loss.

- `'normgini'` = Normalized Gini index.

- `'cont.boyce'` = Continuous Boyce index. Ratio of predicted against
  expected frequency calculated over a moving window:
  \$\$\frac{P\_{i}}{E\_{i}}\$\$ where \\P\_{i} =
  \frac{p\_{i}}{\sum\_{j=1}^{b} p\_{j}}\\ and \\E\_{i} =
  \frac{a\_{i}}{\sum\_{j=1}^{b} a\_{j}}\\.

**Discrete:**

- `'n'` = Number of observations.

- `'auc'` = Area under the curve (AUC), i.e. the integral of a function
  relating the true positive rate against the false positive rate.

- `'overall.accuracy'` = Overall Accuracy: \$\$Accuracy = \frac{TP +
  TN}{N}\$\$

- `'true.presence.ratio'` = True presence ratio or Jaccard index: \$\$J
  = \frac{TP}{TP + TN + FP + FN}\$\$

- `'precision'` = Precision, positive detection rate: \$\$Precision =
  \frac{TP}{TP + FP}\$\$

- `'sensitivity'` = Sensitivity, ratio of true positives against all
  positives: \$\$Sensitivity = \frac{TP}{TP + FN}\$\$

- `'specificity'` = Specificity, ratio of true negatives against all
  negatives: \$\$Specificity = \frac{TN}{TN + FP}\$\$

- `'tss'` = True Skill Statistic: \$\$TSS = Sensitivity + Specificity -
  1\$\$

- `'f1'` = F1 Score or positive predictive value: \$\$F1 = \frac{2 \cdot
  TP}{2 \cdot TP + FP + FN}\$\$

- `'logloss'` = Log loss.

- `'expected.accuracy'` = Expected Accuracy: \$\$EA = \frac{(TP +
  FP)(TP + FN)}{N^2} + \frac{(TN + FN)(TN + FP)}{N^2}\$\$

- `'kappa'` = Cohen's Kappa: \$\$\kappa = \frac{2(TP \cdot TN - FN \cdot
  FP)}{(TP + FP)(FP + TN) + (TP + FN)(FN + TN)}\$\$

- `'brier.score'` = Brier score: \$\$BS = \frac{\sum\_{i=1}^{N}
  (y\_{i} - x\_{i})^{2}}{N}\$\$ where \\y\_{i}\\ is the predicted and
  \\x\_{i}\\ the observed value.

## Note

If you use the Boyce Index, please cite the original Hirzel et al.
(2006) paper.

## References

- Allouche O., Tsoar A., Kadmon R., (2006). Assessing the accuracy of
  species distribution models: prevalence, kappa and the true skill
  statistic (TSS). Journal of Applied Ecology, 43(6), 1223–1232.

- Liu, C., White, M., Newell, G., 2013. Selecting thresholds for the
  prediction of species occurrence with presence-only data. J. Biogeogr.
  40, 778–789. https://doi.org/10.1111/jbi.12058

- Hirzel, A. H., Le Lay, G., Helfer, V., Randin, C., & Guisan, A.
  (2006). Evaluating the ability of habitat suitability models to
  predict species presences. Ecological modelling, 199(2), 142-152.

## Examples

``` r
if (FALSE) { # \dontrun{
 # Assuming that mod is a distribution object and has a thresholded layer
 mod <- threshold(mod, method = "TSS")
 validate(mod, method = "discrete")
 } # }
```

# Summarises a trained model or predictor object

This helper function summarizes a given object, including
[DistributionModel](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md),
[PredictorDataset](https://iiasa.github.io/ibis.iSDM/reference/PredictorDataset-class.md)
or
[PriorList](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
objects and others. This can be a helpful way to summarize what is
contained within and the values of specified models or objects.

When unsure, it is usually a good strategy to run summary on any object.

## Usage

``` r
# S3 method for class 'distribution'
summary(object, ...)

# S3 method for class 'DistributionModel'
summary(object, ...)

# S3 method for class 'PredictorDataset'
summary(object, ...)

# S3 method for class 'BiodiversityScenario'
summary(object, ...)

# S3 method for class 'PriorList'
summary(object, ...)

# S3 method for class 'Settings'
summary(object, ...)
```

## Arguments

- object:

  Any prepared object.

- ...:

  not used.

## See also

[`base::summary()`](https://rdrr.io/r/base/summary.html).

## Examples

``` r
if (FALSE) { # \dontrun{
# Example with a trained model
x <- distribution(background) |>
        # Presence-absence data
        add_biodiversity_poipa(surveydata) |>
        # Add predictors and scale them
        add_predictors(env = predictors) |>
        # Use glmnet and lasso regression for estimation
        engine_glmnet(alpha = 1)
 # Train the model
 mod <- train(x)
 summary(mod)

 # Example with a prior object
 p1 <- BREGPrior(variable = "forest", hyper = 2, ip = NULL)
 p2 <- BREGPrior(variable = "cropland", hyper = NULL, ip = 1)
 pp <- priors(p1,p2)
 summary(pp)
} # }

```

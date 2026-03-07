# Create a posterior prediction from a rstanfit object

This function does simulates from the posterior of a created stan model,
therefore providing a fast and efficient way to project coefficients
obtained from Bayesian models to new/novel contexts.

## Usage

``` r
posterior_predict_stanfit(
  obj,
  form,
  newdata,
  type = "predictor",
  family = NULL,
  offset = NULL,
  draws = NULL
)
```

## Arguments

- obj:

  A `"stanfit"` object (as used by rstan).

- form:

  A [`formula`](https://rdrr.io/r/stats/formula.html) object created for
  the
  [DistributionModel](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md).

- newdata:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) with new data
  to be used for prediction.

- type:

  A [`character`](https://rdrr.io/r/base/character.html) of whether the
  linear `predictor` or the `response` is to be summarized.

- family:

  A [`character`](https://rdrr.io/r/base/character.html) giving the
  family for simulating linear response values (Default: `NULL`)

- offset:

  A [vector](https://rdrr.io/r/base/vector.html) with an optionally
  specified offset.

- draws:

  [numeric](https://rdrr.io/r/base/numeric.html) indicating whether a
  specific number of draws should be taken.

## References

- <https://medium.com/@alex.pavlakis/making-predictions-from-stan-models-in-r-3e349dfac1ed>.

- The brms R-package.

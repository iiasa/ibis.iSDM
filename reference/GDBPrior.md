# Monotonic constrained priors for boosted regressions

Monotonic constrains for gradient descent boosting models do not work in
the same way as other priors where a specific coefficient or magnitude
of importance is specified. Rather monotonic constraints **enforce** a
specific directionality of regression coefficients so that for instance
a coefficient has to be positive or negative.

**Important:** Specifying a monotonic constrain for the
[engine_gdb](https://iiasa.github.io/ibis.iSDM/reference/engine_gdb.md)
does not guarantee that the variable is retained in the model as it can
still be regularized out.

## Usage

``` r
GDBPrior(variable, hyper = "increasing", ...)

# S4 method for class 'character'
GDBPrior(variable, hyper = "increasing", ...)
```

## Arguments

- variable:

  A [`character`](https://rdrr.io/r/base/character.html) matched against
  existing predictors variables.

- hyper:

  A [`character`](https://rdrr.io/r/base/character.html) object
  describing the type of constrain. Available options are
  `'increasing'`, `'decreasing'`, `'convex'`, `'concave'`, `'positive'`,
  `'negative'` or `'none'`.

- ...:

  Variables passed on to prior object.

## Note

Similar priors can also be defined for the
[`engine_xgboost`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md)
via
[`XGBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md).

## References

- Hofner, B., Müller, J., & Hothorn, T. (2011). Monotonicity‐constrained
  species distribution models. Ecology, 92(10), 1895-1901.

## See also

[`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md),
[`XGBPrior`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md)

Other prior:
[`BARTPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPrior.md),
[`BARTPriors()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPriors.md),
[`BREGPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BREGPrior.md),
[`BREGPriors()`](https://iiasa.github.io/ibis.iSDM/reference/BREGPriors.md),
[`GDBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/GDBPriors.md),
[`GLMNETPrior()`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPrior.md),
[`GLMNETPriors()`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPriors.md),
[`INLAPrior()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPrior.md),
[`INLAPriors()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPriors.md),
[`STANPrior()`](https://iiasa.github.io/ibis.iSDM/reference/STANPrior.md),
[`STANPriors()`](https://iiasa.github.io/ibis.iSDM/reference/STANPriors.md),
[`XGBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md),
[`XGBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPriors.md),
[`add_priors()`](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md),
[`get_priors()`](https://iiasa.github.io/ibis.iSDM/reference/get_priors.md),
[`priors()`](https://iiasa.github.io/ibis.iSDM/reference/priors.md),
[`rm_priors()`](https://iiasa.github.io/ibis.iSDM/reference/rm_priors.md)

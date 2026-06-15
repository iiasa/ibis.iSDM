# Creates a new PriorList object

A
[`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
object is essentially a list that contains individual
[`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
objects. In order to use priors for any of the engines, the respective
[`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
has to be identified (e.g.
[`INLAPrior`](https://iiasa.github.io/ibis.iSDM/reference/INLAPrior.md))
and embedded in a
[`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
object. Afterwards these objects can then be added to a
[distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object with the
[add_priors](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md)
function.

A
[`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
object is essentially a list that contains individual
[`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
objects. In order to use priors for any of the engines, the respective
[`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
has to be identified (e.g.
[`INLAPrior`](https://iiasa.github.io/ibis.iSDM/reference/INLAPrior.md))
and embedded in a
[`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
object. Afterwards these objects can then be added to a
[distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object with the
[add_priors](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md)
function.

## Usage

``` r
priors(x, ...)

# S4 method for class 'ANY'
priors(x, ...)

priors(x, ...)

# S4 method for class 'ANY'
priors(x, ...)
```

## Arguments

- x:

  A
  [`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
  object added to the list.

- ...:

  One or multiple additional
  [`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
  object added to the list.

## Value

A
[`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
object.

A
[`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
object.

## See also

[`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md),
[`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)

Other prior:
[`BARTPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPrior.md),
[`BARTPriors()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPriors.md),
[`BREGPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BREGPrior.md),
[`BREGPriors()`](https://iiasa.github.io/ibis.iSDM/reference/BREGPriors.md),
[`GDBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/GDBPrior.md),
[`GDBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/GDBPriors.md),
[`GLMNETPrior()`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPrior.md),
[`GLMNETPriors()`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPriors.md),
[`INLAPrior()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPrior.md),
[`INLAPriors()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPriors.md),
[`STANPrior()`](https://iiasa.github.io/ibis.iSDM/reference/STANPrior.md),
[`STANPriors()`](https://iiasa.github.io/ibis.iSDM/reference/STANPriors.md),
[`XGBInteractionPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBInteractionPrior.md),
[`XGBInteractionPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBInteractionPriors.md),
[`XGBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md),
[`XGBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPriors.md),
[`add_priors()`](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md),
[`get_priors()`](https://iiasa.github.io/ibis.iSDM/reference/get_priors.md),
[`rm_priors()`](https://iiasa.github.io/ibis.iSDM/reference/rm_priors.md)

## Examples

``` r
p1 <- GDBPrior(variable = "Forest", hyper = "positive")
p2 <- GDBPrior(variable = "Urban", hyper = "decreasing")
priors(p1, p2)
#> Set priors: 2

if (FALSE) { # \dontrun{
p1 <- INLAPrior(variable = "Forest",type = "normal", hyper = c(1,1e4))
p2 <- INLAPrior(variable = "Urban",type = "normal", hyper = c(0,1e-2))
priors(p1, p2)
} # }
```

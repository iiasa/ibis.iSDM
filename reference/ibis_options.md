# Print ibis options

There are a number of hidden options that can be specified for
ibis.iSDM. Currently supported are:

- `'ibis.runparallel'` :
  [`logical`](https://rdrr.io/r/base/logical.html) value on whether
  processing should be run in parallel.

- `'ibis.nthread'` : [`numeric`](https://rdrr.io/r/base/numeric.html)
  value on how many cores should be used by default.

- `'ibis.setupmessages'` :
  [`logical`](https://rdrr.io/r/base/logical.html) value indicating
  whether message during object creation should be shown (Default:
  `NULL`).

- `'ibis.engines'` : Returns a
  [`vector`](https://rdrr.io/r/base/vector.html) with all valid engines.

- `'ibis.use_future'` : [`logical`](https://rdrr.io/r/base/logical.html)
  on whether the future package should be used for parallel computing.

## Usage

``` r
ibis_options()
```

## Value

The output of `getOptions` for all ibis related variables.

## Examples

``` r
ibis_options()
#> $ibis.cleannames
#> [1] TRUE
#> 
#> $ibis.corPred
#> [1] 0.7
#> 
#> $ibis.dependencies
#>  [1] "scales"        "biscale"       "modEvA"        "dplyr"        
#>  [5] "geodist"       "geosphere"     "progress"      "glmnet"       
#>  [9] "glmnetUtils"   "xgboost"       "BoomSpikeSlab" "INLA"         
#> [13] "inlabru"       "gnlm"          "cubelyr"       "matrixStats"  
#> [17] "Boruta"        "abess"         "gdalUtilities" "dbarts"       
#> [21] "mboost"        "rstan"         "cmdstanr"      "biscale"      
#> [25] "poems"         "BiocManager"  
#> 
#> $ibis.engines
#>  [1] "GDB-Model"     "BART-Model"    "INLABRU-Model" "BREG-Model"   
#>  [5] "GLMNET-Model"  "GLM-Model"     "SCAMPR-Model"  "INLA-Model"   
#>  [9] "STAN-Model"    "XGBOOST-Model"
#> 
#> $ibis.futurestrategy
#> [1] "sequential"
#> 
#> $ibis.nthread
#> [1] 3
#> 
#> $ibis.priors
#> [1] "INLAPrior"   "BARTPrior"   "GDBPrior"    "GLMNETPrior" "XGBPrior"   
#> [6] "BREGPrior"   "STANPrior"  
#> 
#> $ibis.pseudoabsence
#> Background Settings: 5 parameters
#> 
#> $ibis.runparallel
#> [1] FALSE
#> 
#> $ibis.seed
#> [1] 12830
#> 
#> $ibis.setupmessages
#> [1] TRUE
#> 
```

# Biodiversity Distribution master class

Base [`R6`](https://r6.r-lib.org/reference/R6Class.html) class for any
biodiversity distribution objects. Serves as container that supplies
data and functions to other
[`R6`](https://r6.r-lib.org/reference/R6Class.html) classes. Generally
stores all objects and parameters added to a model.

## Value

An [`R6::R6Class`](https://r6.r-lib.org/reference/R6Class.html)
generator object.

## Details

Run [`names()`](https://rdrr.io/r/base/names.html) on a
[`distribution`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object to show all available functions.

## Note

Not implemented yet.

Not implemented yet.

## See also

[`add_limits_extrapolation()`](https://iiasa.github.io/ibis.iSDM/reference/add_limits_extrapolation.md)

[`add_latent_spatial()`](https://iiasa.github.io/ibis.iSDM/reference/add_latent_spatial.md)

[`add_priors()`](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md)

[`add_biodiversity_poipa()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_poipa.md),
[`add_biodiversity_poipo()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_poipo.md),
[`add_biodiversity_polpa()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_polpa.md),
[`add_biodiversity_polpo()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_polpo.md)

[`add_predictors()`](https://iiasa.github.io/ibis.iSDM/reference/add_predictors.md)

[`add_offset()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset.md)

## Public fields

- `background`:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object
  delineating the modelling extent.

- `limits`:

  An optional [`sf`](https://r-spatial.github.io/sf/reference/sf.html)
  object on potential extrapolation limits

- `biodiversity`:

  A
  [`BiodiversityDatasetCollection`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDatasetCollection-class.md)
  object.

- `predictors`:

  A
  [`PredictorDataset`](https://iiasa.github.io/ibis.iSDM/reference/PredictorDataset-class.md)
  object.

- `priors`:

  An optional
  [`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
  object.

- `control`:

  An optional Control object such as for biases.

- `latentfactors`:

  A [`character`](https://rdrr.io/r/base/character.html) on whether
  latentfactors are used.

- `offset`:

  A [`character`](https://rdrr.io/r/base/character.html) on whether
  methods are used.

- `log`:

  An optional
  [`Log`](https://iiasa.github.io/ibis.iSDM/reference/Log-class.md)
  object.

- `engine`:

  A
  [`Engine`](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md)
  object.

## Methods

### Public methods

- [`BiodiversityDistribution$new()`](#method-BiodiversityDistribution-initialize)

- [`BiodiversityDistribution$print()`](#method-BiodiversityDistribution-print)

- [`BiodiversityDistribution$show()`](#method-BiodiversityDistribution-show)

- [`BiodiversityDistribution$name()`](#method-BiodiversityDistribution-name)

- [`BiodiversityDistribution$show_background_info()`](#method-BiodiversityDistribution-show_background_info)

- [`BiodiversityDistribution$set_limits()`](#method-BiodiversityDistribution-set_limits)

- [`BiodiversityDistribution$get_limits()`](#method-BiodiversityDistribution-get_limits)

- [`BiodiversityDistribution$rm_limits()`](#method-BiodiversityDistribution-rm_limits)

- [`BiodiversityDistribution$get_predictor_names()`](#method-BiodiversityDistribution-get_predictor_names)

- [`BiodiversityDistribution$set_latent()`](#method-BiodiversityDistribution-set_latent)

- [`BiodiversityDistribution$get_latent()`](#method-BiodiversityDistribution-get_latent)

- [`BiodiversityDistribution$rm_latent()`](#method-BiodiversityDistribution-rm_latent)

- [`BiodiversityDistribution$get_priors()`](#method-BiodiversityDistribution-get_priors)

- [`BiodiversityDistribution$set_priors()`](#method-BiodiversityDistribution-set_priors)

- [`BiodiversityDistribution$set_biodiversity()`](#method-BiodiversityDistribution-set_biodiversity)

- [`BiodiversityDistribution$set_predictors()`](#method-BiodiversityDistribution-set_predictors)

- [`BiodiversityDistribution$set_engine()`](#method-BiodiversityDistribution-set_engine)

- [`BiodiversityDistribution$get_engine()`](#method-BiodiversityDistribution-get_engine)

- [`BiodiversityDistribution$rm_engine()`](#method-BiodiversityDistribution-rm_engine)

- [`BiodiversityDistribution$get_prior_variables()`](#method-BiodiversityDistribution-get_prior_variables)

- [`BiodiversityDistribution$set_offset()`](#method-BiodiversityDistribution-set_offset)

- [`BiodiversityDistribution$get_offset()`](#method-BiodiversityDistribution-get_offset)

- [`BiodiversityDistribution$rm_offset()`](#method-BiodiversityDistribution-rm_offset)

- [`BiodiversityDistribution$plot_offsets()`](#method-BiodiversityDistribution-plot_offsets)

- [`BiodiversityDistribution$get_offset_type()`](#method-BiodiversityDistribution-get_offset_type)

- [`BiodiversityDistribution$set_control()`](#method-BiodiversityDistribution-set_control)

- [`BiodiversityDistribution$get_control()`](#method-BiodiversityDistribution-get_control)

- [`BiodiversityDistribution$rm_control()`](#method-BiodiversityDistribution-rm_control)

- [`BiodiversityDistribution$plot_bias()`](#method-BiodiversityDistribution-plot_bias)

- [`BiodiversityDistribution$get_log()`](#method-BiodiversityDistribution-get_log)

- [`BiodiversityDistribution$set_log()`](#method-BiodiversityDistribution-set_log)

- [`BiodiversityDistribution$get_extent()`](#method-BiodiversityDistribution-get_extent)

- [`BiodiversityDistribution$get_projection()`](#method-BiodiversityDistribution-get_projection)

- [`BiodiversityDistribution$get_resolution()`](#method-BiodiversityDistribution-get_resolution)

- [`BiodiversityDistribution$rm_predictors()`](#method-BiodiversityDistribution-rm_predictors)

- [`BiodiversityDistribution$rm_priors()`](#method-BiodiversityDistribution-rm_priors)

- [`BiodiversityDistribution$show_biodiversity_length()`](#method-BiodiversityDistribution-show_biodiversity_length)

- [`BiodiversityDistribution$show_biodiversity_equations()`](#method-BiodiversityDistribution-show_biodiversity_equations)

- [`BiodiversityDistribution$get_biodiversity_equations()`](#method-BiodiversityDistribution-get_biodiversity_equations)

- [`BiodiversityDistribution$get_biodiversity_types()`](#method-BiodiversityDistribution-get_biodiversity_types)

- [`BiodiversityDistribution$get_biodiversity_ids()`](#method-BiodiversityDistribution-get_biodiversity_ids)

- [`BiodiversityDistribution$get_biodiversity_names()`](#method-BiodiversityDistribution-get_biodiversity_names)

- [`BiodiversityDistribution$plot()`](#method-BiodiversityDistribution-plot)

- [`BiodiversityDistribution$summary()`](#method-BiodiversityDistribution-summary)

- [`BiodiversityDistribution$clone()`](#method-BiodiversityDistribution-clone)

------------------------------------------------------------------------

### `BiodiversityDistribution$new()`

Initializes the object and creates an BiodiversityDataset by default.

#### Usage

    BiodiversityDistribution$new(background, limits, biodiversity, ...)

#### Arguments

- `background`:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object
  delineating the modelling extent.

- `limits`:

  An optional [`sf`](https://r-spatial.github.io/sf/reference/sf.html)
  object on potential extrapolation limits

- `biodiversity`:

  A
  [`BiodiversityDatasetCollection`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDatasetCollection-class.md)
  object.

- `...`:

  Any other objects

#### Returns

NULL

------------------------------------------------------------------------

### `BiodiversityDistribution$print()`

Looks for and returns the properties of all contained objects.

#### Usage

    BiodiversityDistribution$print()

#### Returns

A message on screen

------------------------------------------------------------------------

### `BiodiversityDistribution$show()`

An alias for print

#### Usage

    BiodiversityDistribution$show()

#### Returns

A message on screen

------------------------------------------------------------------------

### `BiodiversityDistribution$name()`

Returns self-describing name

#### Usage

    BiodiversityDistribution$name()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the name

------------------------------------------------------------------------

### `BiodiversityDistribution$show_background_info()`

Summarizes extent and projection from set background

#### Usage

    BiodiversityDistribution$show_background_info()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the name

------------------------------------------------------------------------

### `BiodiversityDistribution$set_limits()`

Specify new limits to the background

#### Usage

    BiodiversityDistribution$set_limits(x)

#### Arguments

- `x`:

  A [`list`](https://rdrr.io/r/base/list.html) object with method and
  limit type.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$get_limits()`

Get provided limits if set or a waiver

#### Usage

    BiodiversityDistribution$get_limits()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) or waiver.

------------------------------------------------------------------------

### `BiodiversityDistribution$rm_limits()`

Remove limits if set.

#### Usage

    BiodiversityDistribution$rm_limits()

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$get_predictor_names()`

Function for querying predictor names if existing

#### Usage

    BiodiversityDistribution$get_predictor_names()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) vector.

------------------------------------------------------------------------

### `BiodiversityDistribution$set_latent()`

Adding latent factors to the object.

#### Usage

    BiodiversityDistribution$set_latent(type, method = NULL, separate_spde = FALSE)

#### Arguments

- `type`:

  A [`character`](https://rdrr.io/r/base/character.html) with the given
  type.

- `method`:

  A [`character`](https://rdrr.io/r/base/character.html) with a method.

- `separate_spde`:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag whether
  duplicate of SPDE effects are to be created.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$get_latent()`

Get latent factors if found in object.

#### Usage

    BiodiversityDistribution$get_latent()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with those
objects.

------------------------------------------------------------------------

### `BiodiversityDistribution$rm_latent()`

Remove latent factors if found in object.

#### Usage

    BiodiversityDistribution$rm_latent()

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$get_priors()`

Get prior object if found in object.

#### Usage

    BiodiversityDistribution$get_priors()

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$set_priors()`

Specify new prior object. Overwrites existing ones

#### Usage

    BiodiversityDistribution$set_priors(x)

#### Arguments

- `x`:

  A
  [`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
  object.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$set_biodiversity()`

Adds a new biodiversity object to the existing empty collection.

#### Usage

    BiodiversityDistribution$set_biodiversity(id, p)

#### Arguments

- `id`:

  A [`character`](https://rdrr.io/r/base/character.html) or id defining
  this object.

- `p`:

  A
  [`BiodiversityDataset`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDataset-class.md)
  object.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$set_predictors()`

Set a new Predictor object to this object.

#### Usage

    BiodiversityDistribution$set_predictors(x)

#### Arguments

- `x`:

  A
  [`PredictorDataset`](https://iiasa.github.io/ibis.iSDM/reference/PredictorDataset-class.md)
  with predictors for this object.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$set_engine()`

Set a new Engine object to this object.

#### Usage

    BiodiversityDistribution$set_engine(x)

#### Arguments

- `x`:

  A
  [`Engine`](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md)
  for this object.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$get_engine()`

Gets the name of the current engine if set.

#### Usage

    BiodiversityDistribution$get_engine()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the engine
name

------------------------------------------------------------------------

### `BiodiversityDistribution$rm_engine()`

Removes the current engine if set.

#### Usage

    BiodiversityDistribution$rm_engine()

#### Returns

This object

------------------------------------------------------------------------

### `BiodiversityDistribution$get_prior_variables()`

Get prior variables

#### Usage

    BiodiversityDistribution$get_prior_variables()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the variable
names for which priors have been added.

------------------------------------------------------------------------

### `BiodiversityDistribution$set_offset()`

Specify new offsets.

#### Usage

    BiodiversityDistribution$set_offset(x)

#### Arguments

- `x`:

  A new
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object to be used as offset.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$get_offset()`

Get offset (print name)

#### Usage

    BiodiversityDistribution$get_offset()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with all the
offsets in here.

------------------------------------------------------------------------

### `BiodiversityDistribution$rm_offset()`

Remove offsets if found.

#### Usage

    BiodiversityDistribution$rm_offset(what = NULL)

#### Arguments

- `what`:

  Optional [`character`](https://rdrr.io/r/base/character.html) of
  specific offsets to remove.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$plot_offsets()`

Plot offset if found.

#### Usage

    BiodiversityDistribution$plot_offsets()

#### Returns

A graphical element.

------------------------------------------------------------------------

### `BiodiversityDistribution$get_offset_type()`

Get offset parameters if found

#### Usage

    BiodiversityDistribution$get_offset_type()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) with the offset parameters
if found.

------------------------------------------------------------------------

### `BiodiversityDistribution$set_control()`

Set new bias control

#### Usage

    BiodiversityDistribution$set_control(type = "bias", x, method, value)

#### Arguments

- `type`:

  A [`character`](https://rdrr.io/r/base/character.html) with the type
  of control object.

- `x`:

  A new bias control object. Expecting a
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object.

- `method`:

  The method used to create the object.

- `value`:

  A vector of supplied values. For `"bias"` those should be
  [`numeric`](https://rdrr.io/r/base/numeric.html).

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$get_control()`

Get bias control (print name)

#### Usage

    BiodiversityDistribution$get_control(type = "bias")

#### Arguments

- `type`:

  A [`character`](https://rdrr.io/r/base/character.html) with the type
  of control object.

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the control
object if found.

------------------------------------------------------------------------

### `BiodiversityDistribution$rm_control()`

Remove bias controls if found.

#### Usage

    BiodiversityDistribution$rm_control(type)

#### Arguments

- `type`:

  A [`character`](https://rdrr.io/r/base/character.html) with the type
  of control object.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$plot_bias()`

Plot bias variable if set.

#### Usage

    BiodiversityDistribution$plot_bias()

#### Returns

A graphical element.

------------------------------------------------------------------------

### `BiodiversityDistribution$get_log()`

Returns the output filename of the current log object if set.

#### Usage

    BiodiversityDistribution$get_log()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) where the output
is returned.

------------------------------------------------------------------------

### `BiodiversityDistribution$set_log()`

Set a new log object

#### Usage

    BiodiversityDistribution$set_log(x)

#### Arguments

- `x`:

  A [`Log`](https://iiasa.github.io/ibis.iSDM/reference/Log-class.md)
  object.

#### Returns

This object

------------------------------------------------------------------------

### `BiodiversityDistribution$get_extent()`

Get extent

#### Usage

    BiodiversityDistribution$get_extent()

#### Returns

Background extent or NULL.

------------------------------------------------------------------------

### `BiodiversityDistribution$get_projection()`

Get projection from the background in crs format.

#### Usage

    BiodiversityDistribution$get_projection()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) of the projection

------------------------------------------------------------------------

### `BiodiversityDistribution$get_resolution()`

Return resolution of the background object.

#### Usage

    BiodiversityDistribution$get_resolution()

#### Returns

A [`vector`](https://rdrr.io/r/base/vector.html) with the resolution.

------------------------------------------------------------------------

### `BiodiversityDistribution$rm_predictors()`

Remove predictiors. Either all of them or specific ones.

#### Usage

    BiodiversityDistribution$rm_predictors(names)

#### Arguments

- `names`:

  A [`character`](https://rdrr.io/r/base/character.html) with the
  predictors to be removed.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$rm_priors()`

Remove priors. Either all of them or specific ones.

#### Usage

    BiodiversityDistribution$rm_priors(names = NULL)

#### Arguments

- `names`:

  A [`character`](https://rdrr.io/r/base/character.html) with the priors
  to be removed.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityDistribution$show_biodiversity_length()`

Show number of biodiversity records

#### Usage

    BiodiversityDistribution$show_biodiversity_length()

#### Returns

A [`numeric`](https://rdrr.io/r/base/numeric.html) with sum of
biodiversity records

------------------------------------------------------------------------

### `BiodiversityDistribution$show_biodiversity_equations()`

Show Equations of biodiversity records

#### Usage

    BiodiversityDistribution$show_biodiversity_equations()

#### Returns

A message on screen.

------------------------------------------------------------------------

### `BiodiversityDistribution$get_biodiversity_equations()`

Get equations of biodiversity records

#### Usage

    BiodiversityDistribution$get_biodiversity_equations()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) vector.

------------------------------------------------------------------------

### `BiodiversityDistribution$get_biodiversity_types()`

Query all biodiversity types in this object

#### Usage

    BiodiversityDistribution$get_biodiversity_types()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) vector.

------------------------------------------------------------------------

### `BiodiversityDistribution$get_biodiversity_ids()`

Return all biodiversity dataset ids in the object

#### Usage

    BiodiversityDistribution$get_biodiversity_ids()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) for the ids in the
biodiversity datasets

------------------------------------------------------------------------

### `BiodiversityDistribution$get_biodiversity_names()`

Return all the [`character`](https://rdrr.io/r/base/character.html)
names of all biodiversity datasets

#### Usage

    BiodiversityDistribution$get_biodiversity_names()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) with the names in the
biodiversity datasets

------------------------------------------------------------------------

### `BiodiversityDistribution$plot()`

Plots the content of this class.

#### Usage

    BiodiversityDistribution$plot()

#### Returns

A message.

------------------------------------------------------------------------

### `BiodiversityDistribution$summary()`

Summary function for this object.

#### Usage

    BiodiversityDistribution$summary()

#### Returns

A message.

------------------------------------------------------------------------

### `BiodiversityDistribution$clone()`

The objects of this class are cloneable with this method.

#### Usage

    BiodiversityDistribution$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# Query available functions and entries
background <- terra::rast(system.file('extdata/europegrid_50km.tif',
package='ibis.iSDM',mustWork = TRUE))
# Define model
x <- distribution(background)
#> [Setup] 2026-07-30 07:06:18.110693 | Creating distribution object...
names(x)
#>  [1] ".__enclos_env__"             "engine"                     
#>  [3] "log"                         "offset"                     
#>  [5] "latentfactors"               "control"                    
#>  [7] "priors"                      "predictors"                 
#>  [9] "biodiversity"                "limits"                     
#> [11] "background"                  "clone"                      
#> [13] "summary"                     "plot"                       
#> [15] "get_biodiversity_names"      "get_biodiversity_ids"       
#> [17] "get_biodiversity_types"      "get_biodiversity_equations" 
#> [19] "show_biodiversity_equations" "show_biodiversity_length"   
#> [21] "rm_priors"                   "rm_predictors"              
#> [23] "get_resolution"              "get_projection"             
#> [25] "get_extent"                  "set_log"                    
#> [27] "get_log"                     "plot_bias"                  
#> [29] "rm_control"                  "get_control"                
#> [31] "set_control"                 "get_offset_type"            
#> [33] "plot_offsets"                "rm_offset"                  
#> [35] "get_offset"                  "set_offset"                 
#> [37] "get_prior_variables"         "rm_engine"                  
#> [39] "get_engine"                  "set_engine"                 
#> [41] "set_predictors"              "set_biodiversity"           
#> [43] "set_priors"                  "get_priors"                 
#> [45] "rm_latent"                   "get_latent"                 
#> [47] "set_latent"                  "get_predictor_names"        
#> [49] "rm_limits"                   "get_limits"                 
#> [51] "set_limits"                  "show_background_info"       
#> [53] "name"                        "show"                       
#> [55] "print"                       "initialize"                 
```

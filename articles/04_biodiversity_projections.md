# Creating biodiversity projections

Besides inferring environmental responses and mapping the contemporary
distribution of any species (see [Train simple
Model](https://iiasa.github.io/ibis.iSDM/articles/01_train_simple_model.md)),
species distribution models (SDMs) are commonly used to make predictions
on how biodiversity is expected to change under future conditions.
Usually this is being done by first training a model on present
observations and then projecting the resulting $`\beta`$ coefficients to
another set of environmental predictors from another time period or
spatial region.

The `ibis.iSDM` R-package provides direct support for creating such
projections. This requires a previously fitted \[`DistributionModel`\]
and a new set of covariates which match in names the covariates on which
the fitted model was trained. The key functions here are
[`scenario()`](https://iiasa.github.io/ibis.iSDM/reference/scenario.md)
and
[`project()`](https://iiasa.github.io/ibis.iSDM/reference/project.md)
which are specific to such projections and result in the creation of a
\[`BiodiversityScenario`\] object. For other functions the same syntax
as for previously trained models applies, e.g. new predictor can be
added via
[`add_predictors()`](https://iiasa.github.io/ibis.iSDM/reference/add_predictors.md)
or thresholds applied via
[`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md).
In addition the `ibis.iSDM` R-package allows the specification of a
number of constraints on the projections, such as for instance dispersal
constraints based on the expected or simulated chance of individuals
dispersing to neighbouring grid cells. Such constraints can be added via
[`add_constraint()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint.md)
(and other constraint functionalities in `add_constraint_*()`).

**Note**: Projecting relative suitability estimates into conditions
other than the observational records for which the SDM was trained comes
with a number of assumptions, most importantly that relationships
between occurrences and environmental conditions are in equilibrium and
similar biases and conditions apply in the future. See [Elith et
al. (2010)](https://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2010.00036.x)
and discussion in [Zurell et
al. (2016)](http://doi.wiley.com/10.1111/gcb.13251) for an introduction
and comparative overview.

## Load relevant packages and testing data

``` r

# Load the packages
library(ibis.iSDM)
library(stars)
library(xgboost)
library(terra)
library(igraph)
library(ggplot2)
library(ncdf4)
library(assertthat)

# Don't print out as many messages
options("ibis.setupmessages" = FALSE)
```

For the purpose of this example we are loading some testing data on
species distributions as well as contemporary and future predictors.
Note the names of predictors used for building a distribution model have
to be consistent with those for creating projections!

``` r

# Background and biodiversity data
background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM'))
virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM'), 'points', quiet = TRUE)

# Note we are loading different predictors than in previous examples
# These are in netcdf4 format, a format specific for storing spatial-temporal data including metadata.
ll <- list.files(system.file("extdata/predictors_presfuture/", package = "ibis.iSDM", mustWork = TRUE), "*.nc",full.names = TRUE)

# From those list of predictors are first loading the current ones as raster data
# We are loading only data from the very first, contemporary time step for model fitting
pred_current <- terra::rast()
for(i in ll) suppressWarnings( pred_current <- c(pred_current, terra::rast(i, lyrs = 1) ) )
names(pred_current) <- tools::file_path_sans_ext( basename(ll) )

# Get future predictors
# These we will load in using the stars package and also ignoring the first time step
pred_future <- stars::read_stars(ll) |> stars:::slice.stars('Time', 2:86)
sf::st_crs(pred_future) <- sf::st_crs(4326) # Set projection
# Rename future predictors to those of current
names(pred_future) <- names(pred_current)

# Plot the test data
plot(pred_current['secdf'],
     col = colorRampPalette(c("grey20", "orange", "lightgreen", "green"))(10),
     main = "Share of secondary vegetation")
```

![](04_biodiversity_projections_files/figure-html/Load%20data-1.png) \##
Train model and create a future projection

We will make use of the data loaded above to (a) first create a species
distribution model for contemporary conditions and (b) project the
obtained coefficients into the future using future predictors. For
guidance on how distribution models are trained, see other vignettes
([1](https://iiasa.github.io/ibis.iSDM/articles/01_train_simple_model.md)).

``` r

# Train model adding the data loaded above
x <- distribution(background) |> 
  add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed', name = 'Virtual points') |> 
  # Note that we scale the predictors here
  add_predictors(pred_current, transform = 'scale',derivates = 'none') |> 
  engine_glmnet(alpha = 0) 
#> Loaded glmnet 5.0

# Train the model
modf <- train(x, runname = 'Simple PPM', verbose = FALSE)

# Add a threshold to this model by getting 05 percentile of values
modf <- threshold(modf, method = 'percentile', value = 0.05)

# -- #
# Now lets create a scenarios object via scenarios
sc <- scenario(modf) |> 
  # Apply the same variable transformations as above. 
  add_predictors(pred_future, transform = 'scale') |> 
  # Calculate thresholds at each time step. The threshold estimate is taken from the model object.
  threshold()
#> ! State variable of transformation not found?

# This creates a scenario object
sc
#> Spatial-temporal scenario:
#>   Used model: GLMNET-Model
#>  --------- 
#>   Predictors:     bio01, bio12, crops, ... (9 predictors)
#>   Time period:    2016-01-01 -- 2100-01-01 (83.9 years)
#>  --------- 
#>   Threshold:      0.031 (percentile)
#>  --------- 
#>   Scenarios fitted: None
# The object contains its own functions. See the scenarios help file for more information on 
# what is possible with them
names(sc)
#>  [1] "threshold"            "verify"               "summary_beforeafter" 
#>  [4] "summary"              "show"                 "set_simulation"      
#>  [7] "set_predictors"       "set_log"              "set_latent"          
#> [10] "set_data"             "set_constraints"      "scenarios"           
#> [13] "save"                 "rm_predictors"        "rm_limits"           
#> [16] "rm_latent"            "rm_data"              "rm_constraints"      
#> [19] "print"                "predictors"           "plot_threshold"      
#> [22] "plot_scenarios_slope" "plot_relative_change" "plot_migclim"        
#> [25] "plot_animation"       "plot"                 "modelobject"         
#> [28] "modelid"              "mask"                 "log"                 
#> [31] "limits"               "latentfactors"        "initialize"          
#> [34] "get_timeperiod"       "get_thresholdvalue"   "get_threshold"       
#> [37] "get_simulation"       "get_resolution"       "get_projection"      
#> [40] "get_predictors"       "get_predictor_names"  "get_model"           
#> [43] "get_log"              "get_limits"           "get_latent"          
#> [46] "get_data"             "get_constraints"      "get_centroid"        
#> [49] "constraints"          "clone"                "calc_scenarios_slope"
#> [52] "apply_threshold"      ".__enclos_env__"
```

The scenario object can finally be trained via
[`project()`](https://iiasa.github.io/ibis.iSDM/reference/project.md).

``` r

sc.fit1 <- sc |> project()
# Note that an indication of fitted scenarios has been added to the object
sc.fit1
#> Spatial-temporal scenario:
#>   Used model: GLMNET-Model
#>  --------- 
#>   Predictors:     bio01, bio12, crops, ... (9 predictors)
#>   Time period:    2016-01-01 -- 2100-01-01 (83.9 years)
#>  --------- 
#>   Threshold:      0.031 (percentile)
#>  --------- 
#>   Scenarios fitted: Yes
```

## Summarizing and plotting the fitted projections

As with distribution models there are a number of ways how the scenarios
can be visualized and interacted with:

- [`plot()`](https://iiasa.github.io/ibis.iSDM/reference/plot.md) makes
  a visualization of the projections over all time steps (!)
- `plot_relative_change()` calculates the change in suitability area
  between the first and the last timestep and categorizes the result
  accordingly. Note that SDMs as such cannot directly infer colonization
  or extinction, but only gains or losses of suitable habitat!
- `calc_scenarios_slope()` calculates the slope (rate of change) across
  timesteps. Useful for summarizing results
- [`summary()`](https://iiasa.github.io/ibis.iSDM/reference/summary.md)
  creates a summary output of the contained scenarios. If a
  [`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)
  is specified, this function will summarize the amount of area at each
  timestep.
- [`get_data()`](https://iiasa.github.io/ibis.iSDM/reference/get_data.md)
  gets the created scenarios a `stars` object (plus *thresholds* if
  specified).

``` r

# Plot all scenarios. With a large number of predictors this figure will be messy...
plot(sc.fit1) # or sc.fit1$plot()
```

![](04_biodiversity_projections_files/figure-html/Plotting%20and%20summarizing%20the%20created%20projections-1.png)

``` r


# As an alternative, visualize the linear slope per grid cell and across all time steps
o <- sc.fit1$calc_scenarios_slope(plot = TRUE)
```

![](04_biodiversity_projections_files/figure-html/Plotting%20and%20summarizing%20the%20created%20projections-2.png)

``` r


# Another option is to calculate the relative change between start and finish
o <- sc.fit1$plot_relative_change(plot = TRUE)
```

![](04_biodiversity_projections_files/figure-html/Plotting%20and%20summarizing%20the%20created%20projections-3.png)

``` r


# We can also summarize the thresholded data
o <- sc.fit1$summary()
plot(area_km2~band, data = o, type = 'b',
     main = "Suitable habitat across Time",
     ylab = "Amount of area (km2)", xlab = "Time")
```

![](04_biodiversity_projections_files/figure-html/Plotting%20and%20summarizing%20the%20created%20projections-4.png)

``` r


# How does habitat gain and loss change over time?
plot(totchange_gain_km2~band, data = o, type = 'n',
     main = "Habitat gain and loss", ylim = c(-1.5e4, 1.5e4),
     ylab = "Amount of area (km2)", xlab = "Time")
lines(o$totchange_gain_km2~o$band, col = "blue")
lines((o$totchange_loss_km2)~o$band, col = "red")
```

![](04_biodiversity_projections_files/figure-html/Plotting%20and%20summarizing%20the%20created%20projections-5.png)

Finally, scenarios projections can also be saved as specific outputs. As
before, this is enabled via
[`write_output()`](https://iiasa.github.io/ibis.iSDM/reference/write_output.md)
and works just the same for \[`BiodiversityScenario`\] objects, with the
only difference being that the output can be specified as netCDF-4 file.

## Adding constraints to projections

In the simple scenario above we use the naive assumption that any,
depending on the response functions of the fitted distribution model,
any suitable habitat within the background modelling region is
potentially reachable by the species. In reality there might however be
geographic (e.g. islands), environmental and biotic constraints on how
far a species can disperse. These can be specified through the constrain
function
\[[`add_constraint()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint.md)\]
and a variety of constraints is currently available, some of which
depend on other packages.

- [`add_constraint()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint.md)
  Generic wrapper to which a specific ‘method’ can be supplied. See
  documentation for more information on available options and
  parameters.
- [`add_constraint_dispersal()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_dispersal.md)
  To add a dispersal constraint to the projections which is applied
  after each time step. Supports various options with `'sdd_fixed'` for
  fixed dispersal kernels, `'sdd_nexpkernel'` for a negative exponential
  kernel or `'kissmig'` for applying the kissmig framework.
- [`add_constraint_MigClim()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_MigClim.md)
  Use the MigClim R-package to simulate dispersal events between time
  steps. A number of parameters are required here and adding this
  constrain will also overwrite some default plotting capacities (For
  example via `sc$plot_migclim()`). See also the help file and [Engler
  et
  al. (2012)](https://onlinelibrary.wiley.com/doi/10.1111/j.1600-0587.2012.07608.x)
  for more information.
- [`add_constraint_connectivity()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_connectivity.md)
  Add a connectivity constrain to the projection. Currently only hard
  barriers are implemented, but in future additional sub-modules are
  planned to enable more options here.
- [`add_constraint_adaptability()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_adaptability.md)
  Simple constraints on the adaptability of species to novel climatic
  conditions. Currently only simple nichelimits are implemented, which
  ‘cap’ projections in novel environments to the observed ranges of
  contemporary predictors.
- [`add_constraint_boundary()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_boundary.md)
  Specifying a hard boundary constraint on all projections, for example
  by limiting (future) projections to only a certain area such as a
  biome or a contemporary range.

Lastly there are also options to *stabilize* suitability projections via
the
[`project()`](https://iiasa.github.io/ibis.iSDM/reference/project.md)
function. Specifying a stabilization here results in the projections
being smoothed and informed of incremental time steps. This can
particularly help for projections that use variables known to make
sudden, abrupt jumps between time steps (e.g. precipitation anomalies).

``` r

# Adding a simple negative exponential kernel to constrain the predictions
sc.fit2 <- sc |>   
   add_constraint(method = "sdd_nex", value = 1e5) |> 
   # Directly fit the object
   project(stabilize = F)

# Also fit one projection a nichelimit has been added
sc.fit3 <- sc |>   
   add_constraint(method = "sdd_nex", value = 1e5) |>
   add_constraint_adaptability(method = "nichelimit") |> 
   # Directly fit the object
   project(stabilize = F)

# Note how constrains are indicated in the scenario object.
sc.fit3
#> Spatial-temporal scenario:
#>   Used model: GLMNET-Model
#>  --------- 
#>   Predictors:     bio01, bio12, crops, ... (9 predictors)
#>   Time period:    2016-01-01 -- 2100-01-01 (83.9 years)
#>  --------- 
#>   Constraints:    dispersal (sdd_nexpkernel), adaptability (nichelimit)
#>   Threshold:      0.031 (percentile)
#>  --------- 
#>   Scenarios fitted: Yes

# The naive assumption is that there is unlimited dispersal across the whole background
# Note how the projection with dispersal constrain results in a considerable smaller amount of suitable habitat.
sc.fit1$plot(which = 40) # Baseline
```

![](04_biodiversity_projections_files/figure-html/Add%20constraints%20and%20reproject-1.png)

``` r

sc.fit2$plot(which = 40) # With dispersal constrain
```

![](04_biodiversity_projections_files/figure-html/Add%20constraints%20and%20reproject-2.png)

``` r

sc.fit3$plot(which = 40) # With dispersal limit and nichelimitation (within a standard deviation)
```

``` r

# Lets compare the difference in projections compared to the naive one defined earlier. 
o1 <- sc.fit1$summary()
o2 <- sc.fit2$summary()
o3 <- sc.fit3$summary()
arlim <- c(min(o1$area_km2, o2$area_km2, o3$area_km2)-10000,
           max(o1$area_km2, o2$area_km2, o3$area_km2))

plot(area_km2~band, data = o1, type = 'n',
     ylim = arlim,
     main = "Suitable habitat projection",
     ylab = "Amount of area (km2)", xlab = "Time")
lines(o1$area_km2~o1$band, col = "black", lty = 1)
lines(o2$area_km2~o2$band, col = "black", lty = 2)
lines(o3$area_km2~o3$band, col = "black", lty = 3)
legend("bottomleft", 
  legend = c("Unlimited dispersal", "Constrained dispersal",
             "Constrained dispersal and niche limit"), 
  lty = c(1, 2, 3),
  cex = 1.2,
  bty = "n")
```

![](04_biodiversity_projections_files/figure-html/Summarize%20area%20of%20projections-1.png)

``` r


# Lastly it is also possible to directly summarize the state 
# before (usually first year) and end (last year).
sc.fit2$summary_beforeafter()
#> # A tibble: 13 × 5
#>    runname    category                  period        value unit      
#>    <chr>      <chr>                     <chr>         <dbl> <chr>     
#>  1 Simple PPM Current range             2016-01-01  433.    ha        
#>  2 Simple PPM Future range              2100-01-01  329.    ha        
#>  3 Simple PPM Unsuitable                84 years    859.    ha        
#>  4 Simple PPM Loss                      84 years    103.    ha        
#>  5 Simple PPM Gain                      84 years      0     ha        
#>  6 Simple PPM Stable                    84 years    329.    ha        
#>  7 Simple PPM Percent loss              84 years     23.9   %         
#>  8 Simple PPM Percent gain              84 years      0     %         
#>  9 Simple PPM Range change              84 years   -103.    ha        
#> 10 Simple PPM Percent change            84 years    -10.7   %         
#> 11 Simple PPM Sorensen index            84 years      0.875 similarity
#> 12 Simple PPM Centroid distance         84 years    116.    km        
#> 13 Simple PPM Centroid change direction 84 years     32.3   deg
```

Another option for constraining prediction is also by imposing a zonal
limit (for instance climatically defined) on the projections (see
alternatively
[`add_constraint_boundary()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_boundary.md)
above). This has to be done while fitting the SDM for the reference
conditions (see the example with limits
([1](https://iiasa.github.io/ibis.iSDM/articles/01_train_simple_model.md))
) and is considered when doing (future) projections.


<!-- README.md is generated from README.Rmd. Please use this file for any edits-->

# The ibis framework - An **I**ntegrated model for **B**iod**I**versity distribution projection**S**

<!-- <a href='https://github.com/iiasa/rN2000'><img src="man/figures/logo.png" align="right"height=140/></a> --->
<!-- https://shields.io/  For Badges later -->

The ibis.iSDM package provides a series of convenience functions to fit
integrated Species Distribution Models (iSDMs). With integrated models
we generally refer to SDMs that incorporate information from different
biodiversity datasets, external parameters such as priors or offsets
with respect to certain variables and regions. See [Fletcher et
al. (2019)](https://onlinelibrary.wiley.com/doi/abs/10.1002/ecy.2710)
and [Isaac et
al. (2020)](https://linkinghub.elsevier.com/retrieve/pii/S0169534719302551)
for an introduction to iSDMs.

## Installation

Currently only as development version on GitHub.

``` r
install.packages("devtools")
devtools::install_github("IIASA/ibis.iSDM")
```

## Usage

See relevant reference sites and vignettes.

Note that the package in constant development and individual functions
might change. As soon as the package has reached a stable state, the
repository will be made public An official release will be made as soon
as the package is public **In development !!!**

## Package development

(also see [issues](https://github.com/iiasa/ibis.iSDM/issues) and
[projects](https://github.com/iiasa/ibis.iSDM/projects))

-   ✅ Fitting Point Poisson Models (PPMs) using INLA and at least
    another Engine
-   ✅ Basic documentation present that however can be further improved
-   ✅ Various functions to include offsets, ranges, predictor
    transformations
-   ✅ Plotting functions for models and effect are there
-   ✅ Actually support multiple likelihoods depending on included
    biodiversity dataset
-   ✅ Add new class for priors and ensure that they linked in a flexible
    manner
-   🚧 Add a wrapper for PPMs in Stan as new engine for max flexibility
-   🚧 Implement options to add occupancy model for detectability, see
    [Altwegg and
    Nichols](https://onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13090)
    and [Guillera-Arroita](http://doi.wiley.com/10.1111/ecog.02445)
-   🚧 There are reports that iCAR can
    [outperforms](https://arxiv.org/pdf/1204.6087v1.pdf) Matern
    covariance? It likely is also considerably faster, thus could be
    implemented
-   🚧 Create a pkgdown website for the package with tutorials (in
    progress. Looking at it)
-   …

<!-- get_contributors(org = "IIASA", repo = "ibis.iSDM") -->

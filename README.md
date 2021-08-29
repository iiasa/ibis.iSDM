
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
for an introduction.

Isn’t this just another SDM package? Indeed it is, but building our own
modelling infrastructure allows greater flexibility in how models are
created and run. It is less about reinventing the wheel, but rather
bringing strengths from different tools together. *Not exhaustive list
acknowledging other SDM packages in R:*

-   [hSDM](https://github.com/ghislainv/hSDM) -&gt; Bayesian framework
    for mixed models. Fast, but little flexibility with regards to
    weights or offsets
-   [multispeciesPP](https://github.com/wfithian/multispeciesPP) -&gt;
    Not been further developed since years
-   [inlabru](https://github.com/inlabru-org/inlabru) -&gt; Package
    specifically for LGCP with INLA, however little support to other
    engines, likelihoods or data
-   [pointedSDMs](https://github.com/oharar/PointedSDMs) -&gt; INLA
    based SDM package for integrating different datasets.
-   [biomod2](https://github.com/biomodhub/biomod2) -&gt; Popular
    package for ensemble modeling, but no integration, not Bayesian
-   [sdmTMB](https://github.com/pbs-assess/sdmTMB) -&gt; Package for
    fitting spatial-Temporal SDMs
-   [modleR](https://github.com/Model-R/modleR) -&gt; similar as biomod2
    a wrapper to construct ensembles of models
-   [kuenm](https://peerj.com/articles/6281/) -&gt; Another wrapper for
    Maxent

## Installation

Currently only as development version on GitHub.

``` r
install.packages("devtools")
devtools::install_github("IIASA/ibis")
```

## Features and next key milestones

-   ✅ Fitting Point Poisson Models (PPMs) using INLA and at least
    another Engine
-   ✅ Basic documentation present that however can be further improved
-   ✅ Various functions to include offsets, ranges, predictor
    transformations
-   ✅ Plotting functions for models and effect are there
-   🚧 Actually support multiple likelihoods depending on included
    biodiversity dataset
-   🚧 Add a wrapper for PPMs in Stan as new engine for max flexibility
-   🚧 Implement options to add occupancy model for detectability, see
    [Altwegg and
    Nichols](https://onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13090)
    and [Guillera-Arroita](http://doi.wiley.com/10.1111/ecog.02445)
-   🚧 Add new class for priors and ensure that they linked in a flexible
    manner
-   🚧 There are reports that iCAR can
    [outperforms](https://arxiv.org/pdf/1204.6087v1.pdf) Matern
    covariance? It likely is also considerably faster, thus could be
    implemented
-   🚧 Create a pkgdown website for the package with tutorials (see
    below)
-   …

## Development guidelines

-   The ibis.iSDM contains primarily functions for fitting models.
-   Speed and flexibility are key
-   Don’t repeat yourself. Create new functions and if necessary
    classes. Equally try to reuse common names from R, e.g. *plot*,
    *summary*
-   Avoid using additional package dependencies where possible
-   Comment your code!!
-   Use assertions to verify inputs to functions
-   If bored, please write unit tests and ensure they all evaluate
    (CRTL+SHIFT+T)!

*Note that pushing to master or forking is disabled at the moment. Pull
requests require a confirmation.*

## Usage

Note that the package in constant development and individual functions
might change. As soon as the package has reached a stable state, the
repository will be made public An official release will be made as soon
as the package is public **In development !!!**

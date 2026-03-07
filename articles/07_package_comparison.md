# Capabilities compared to other SDM packages

## Comparison with other packages

Species distribution modelling (SDM) approaches have been around for
quite some while and as a result a number of ecological modelling
focused packages have been developed. In general most of them are
customized towards specific purposes and modelling paradigm. So isn’t
this just another SDM package? Indeed it is, but `ibis.iSDM` has a
number of particular features that set it apart from most other SDM
packages:

1.  It focuses particularly on integration as a guiding principle and
    different ways of how heterogeneous sources of evidence can be
    integrated.
2.  It puts a strong focus on biodiversity types, in particular
    Poisson-Process models (PPMs) as the default way of analyzing
    presence-only data.
3.  It follows an object-based and modular programming philosophy,
    taking inspiration from `tidy` programming approaches.
4.  It supports a number of Bayesian SDM approaches and algorithms, a
    field traditionally less represented owing to computational
    constraints.
5.  It is customized to create and modify spatio-temporal scenarios,
    including IIASA integrated land-use assessment model
    [GLOBIOM](https://iiasa.github.io/GLOBIOM/index.html).

Thus overall, the idea package is in part trying to bring some
innovation to the SDM modelling world, while also trying to bring
together strengths from different existing tools.

A Non exhaustive list acknowledging other SDM packages in R and how they
compare to `ibis.iSDM` is provided here:

- [hSDM](https://github.com/ghislainv/hSDM) -\> Bayesian framework for
  hierarchical and mixed models. Fast, but little flexibility with
  regards to weights, offsets and different datatypes.
- [multispeciesPP](https://github.com/wfithian/multispeciesPP) -\>
  Package that allows integrated SDMs, however it has not been further
  developed since a few years and key gaps remain particular with
  regards to different modelling approaches.
- [inlabru](https://github.com/inlabru-org/inlabru) -\> Package
  specifically for Lox-Gaussian-Cox Process (LGCP) models with INLA, now
  integrated also as engine in *ibis.iSDM*
- [pointedSDMs](https://github.com/oharar/PointedSDMs) -\> Another
  wrapper for INLA that allows to integrate different datasets in a SDM.
  Less focus on priors, offsets and scenarios.
- [biomod2](https://github.com/biomodhub/biomod2) -\> Popular package
  for ensemble modelling, but fixed on specific (non-Bayesian) engines
  and data types and with no integration options.
- [sdmTMB](https://github.com/pbs-assess/sdmTMB) -\> Package for fitting
  spatial-Temporal SDMs for very specific biodiversity data.
- [modleR](https://github.com/Model-R/modleR) -\> similar as biomod2 and
  a wrapper to construct ensembles of models.
- [kuenm](https://peerj.com/articles/6281/) -\> Another wrapper for
  Maxent.
- [flexSDM](https://sjevelazco.github.io/flexsdm/) Similar as biomod2 a
  wrapper for SDMs, but coming with several helper functions for data
  preparation and cross-validation.

Besides SDMs there are also new packages available on spatial and
integrated species occupancy models, such as
[spOccupancy](https://www.jeffdoser.com/files/spoccupancy-web/).
Occupancy modelling however requires some very specific biodiversity
data and information to infer detectability of species occurrences.

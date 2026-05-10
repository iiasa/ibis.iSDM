# Comparison of different engines

## Capabilities of included engines

As outlined by [Fletcher et
al. (2019)](https://onlinelibrary.wiley.com/doi/abs/10.1002/ecy.2710),
there are many different forms of integration such as through
\[`ensemble`\] modelling, adding \[`offsets`\], predictors
(e.g. \[[`add_predictor_range()`](https://iiasa.github.io/ibis.iSDM/reference/add_predictor_range.md)\]
) or \[`priors`\] and through full integration of different likelihoods
(See ([Data
integration](https://iiasa.github.io/ibis.iSDM/articles/03_integrate_data.md))
). Not all of these options are available for every engine supported by
the *ibis.iSDM* package and the table below shows the currently
implemented engines and various types of integrations supported by them.

Stating the name and function call of each engine and its supported
model complexity with linear (ln) and non-linear (nl) formulations,
although it should be noted that linear models can approximate
non-linearity by including transformations (as with Maxent,
e.g. hinge/product/quadratic). Not every engine supports the different
types of integration via `ensembles`, `offsets`, `priors`, joint
likelihood estimation and `ensemble` compositing of models using
separate datasets of the same species. When multiple biodiversity
datasets are added to an engine that does not support joint likelihood
estimation, the parameter `method_integration` in
\[[`train()`](https://iiasa.github.io/ibis.iSDM/reference/train.md)\]
determines how the different predictions are integrated. Available
options for integration are via `predictors`, `offsets`, `interactions`,
`priors` or `weights` (see the help file of
\[[`train()`](https://iiasa.github.io/ibis.iSDM/reference/train.md)\]
for more information).

| Name | Complexity | Engine | Offsets | Priors | Weights | Joint likel. |
|----|:--:|:--:|:--:|:--:|:--:|---:|
| Generalized linear model (GLM) | ln | \[[`engine_glm()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glm.md)\] | x |  | x |  |
| Regularized elastic net regression (GLMNET) | ln | \[[`engine_glmnet()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glmnet.md)\] | x | [`GLMNETPrior()`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPrior.md) | x |  |
| Bayesian additive regression trees (BART) | nl | \[[`engine_bart()`](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md)\] | \(x\) | [`BARTPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPrior.md) | x |  |
| Bayesian regularized regression (BREG) | ln | \[[`engine_breg()`](https://iiasa.github.io/ibis.iSDM/reference/engine_breg.md)\] |  | [`BREGPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BREGPrior.md) | x |  |
| Approximate point modelling (SCAMPR) | ln | \[[`engine_scampr()`](https://iiasa.github.io/ibis.iSDM/reference/engine_scampr.md)\] | x |  |  | x |
| Gradient descent boosting (GDB) | ln/nl | \[[`engine_gdb()`](https://iiasa.github.io/ibis.iSDM/reference/engine_gdb.md)\] | x | [`GDBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/GDBPrior.md) | x |  |
| Integrated Nested Laplace approximation (INLA) | ln | \[[`engine_inla()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md)\] | x | [`INLAPrior()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPrior.md) | x | x |
| Integrated Nested Laplace approximation (INLABRU) | ln | \[[`engine_inlabru()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md)\] | x | [`INLAPrior()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPrior.md) | x | x |
| Bayesian regressions (Stan) | ln | \[[`engine_stan()`](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md)\] | x | [`STANPrior()`](https://iiasa.github.io/ibis.iSDM/reference/STANPrior.md) | x | \(x\) |
| eXtreme Gradient Boosting (XGBOOST) | ln/nl | \[[`engine_xgboost()`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md)\] | x | [`XGBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md) | x |  |

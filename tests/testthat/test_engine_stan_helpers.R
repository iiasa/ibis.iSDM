test_that("Stan helpers generate shared-beta code for multiple datasets", {
  build_components <- getFromNamespace("stan_component_table", "ibis.iSDM")
  build_code <- getFromNamespace("stan_model_code", "ibis.iSDM")
  build_data <- getFromNamespace("stan_model_data", "ibis.iSDM")

  settings <- Settings$new()
  settings$set("optim_hyperparam", FALSE)

  model <- list(priors = new_waiver())
  model$biodiversity <- list(
    po1 = list(
      name = "po1", type = "poipo", family = "poisson", link = new_waiver(),
      predictors = data.frame(x = 1:2, y = 1:2, p1 = c(0.1, 0.2), p2 = c(1, 2)),
      predictors_names = c("p1", "p2"),
      observations = data.frame(observed = c(1L, 0L)),
      stan_offset = c(0, 1), stan_weight = c(1, 1), expect = c(1, 1),
      use_intercept = TRUE
    ),
    po2 = list(
      name = "po2", type = "poipo", family = "poisson", link = new_waiver(),
      predictors = data.frame(x = 3:4, y = 3:4, p1 = c(0.3, 0.4), p2 = c(3, 4)),
      predictors_names = c("p1", "p2"),
      observations = data.frame(observed = c(1L, 0L)),
      stan_offset = c(0, 1), stan_weight = c(1, 1), expect = c(1, 1),
      use_intercept = TRUE
    )
  )
  model$stan <- list(
    components = build_components(model),
    feature_names = c("p1", "p2"),
    use_dataset_intercepts = TRUE
  )

  code <- build_code(model, settings)
  dl <- build_data(model, settings)

  expect_match(code, "Integrated multi-dataset species distribution model")
  expect_match(code, "inv_cloglog_response")
  expect_match(code, "poisson_log_lpmf")
  expect_match(code, "vector\\[K\\] beta")
  expect_match(code, "array\\[N\\] int<lower=0> observed")
  expect_false(grepl("\\{\\{", code))
  expect_false(grepl("ibis_inv_cloglog", code))
  expect_false(grepl("vector\\[K\\] b;", code))
  expect_equal(dl$J, 2L)
  expect_equal(dl$J_intercept, 2L)
  expect_equal(dl$K, 2L)
  expect_equal(length(dl$observed), 4L)
})

test_that("Stan model code selects explicit single-dataset templates", {
  build_components <- getFromNamespace("stan_component_table", "ibis.iSDM")
  build_code <- getFromNamespace("stan_model_code", "ibis.iSDM")

  settings <- Settings$new()
  settings$set("optim_hyperparam", FALSE)

  po_model <- list(priors = new_waiver())
  po_model$biodiversity <- list(
    po = list(
      name = "po", type = "poipo", family = "poisson", link = new_waiver(),
      predictors_names = c("p1", "p2"), use_intercept = TRUE
    )
  )
  po_model$stan <- list(
    components = build_components(po_model),
    feature_names = c("p1", "p2"),
    use_dataset_intercepts = TRUE
  )

  pa_model <- list(priors = new_waiver())
  pa_model$biodiversity <- list(
    pa = list(
      name = "pa", type = "poipa", family = "binomial", link = new_waiver(),
      predictors_names = c("p1", "p2"), use_intercept = TRUE
    )
  )
  pa_model$stan <- list(
    components = build_components(pa_model),
    feature_names = c("p1", "p2"),
    use_dataset_intercepts = TRUE
  )

  po_code <- build_code(po_model, settings)
  pa_code <- build_code(pa_model, settings)

  expect_match(po_code, "Single presence-only point process model")
  expect_match(po_code, "poisson_log_lpmf")
  expect_false(grepl("bernoulli", po_code))
  expect_match(po_code, "array\\[N\\] int<lower=0> observed")
  expect_false(grepl("\\{\\{", po_code))

  expect_match(pa_code, "Single presence-absence model")
  expect_match(pa_code, "inv_cloglog_response")
  expect_match(pa_code, "bernoulli_logit_lpmf")
  expect_false(grepl("poisson_log_lpmf", pa_code))
  expect_match(pa_code, "array\\[N\\] int<lower=0> observed")
  expect_false(grepl("\\{\\{", pa_code))
  expect_false(grepl("ibis_inv_cloglog", pa_code))
})

test_that("Stan helpers set PA link defaults by model context", {
  build_components <- getFromNamespace("stan_component_table", "ibis.iSDM")
  resolve_component <- getFromNamespace("stan_component", "ibis.iSDM")
  make_bd <- function(name, type, family, link = new_waiver()) {
    list(name = name, type = type, family = family, link = link, use_intercept = TRUE)
  }

  pa_only <- list(biodiversity = list(pa = make_bd("pa", "poipa", "binomial")))
  pa_only$stan <- list(components = build_components(pa_only))
  expect_equal(pa_only$stan$components$link, "logit")

  mixed <- list(biodiversity = list(
    po = make_bd("po", "poipo", "poisson"),
    pa = make_bd("pa", "poipa", "binomial")
  ))
  mixed$stan <- list(components = build_components(mixed))
  expect_equal(mixed$stan$components$link[2], "cloglog")
  expect_equal(resolve_component(mixed), 1L)
})

test_that("Stan explicit priors are written for shared beta coefficients", {
  build_components <- getFromNamespace("stan_component_table", "ibis.iSDM")
  build_code <- getFromNamespace("stan_model_code", "ibis.iSDM")

  settings <- Settings$new()
  settings$set("optim_hyperparam", FALSE)

  model <- list(priors = priors(STANPrior("p1", "normal", c(1, 0.5))))
  model$biodiversity <- list(
    po = list(
      name = "po", type = "poipo", family = "poisson", link = new_waiver(),
      predictors_names = c("p1", "p2"), use_intercept = TRUE
    )
  )
  model$stan <- list(
    components = build_components(model),
    feature_names = c("p1", "p2"),
    use_dataset_intercepts = TRUE
  )

  code <- build_code(model, settings)
  expect_match(code, "normal_lpdf\\(beta\\[1\\] \\| 1, 0.5\\)")
  expect_match(code, "normal_lpdf\\(beta\\[2\\] \\| 0, 2\\)")
})

test_that("Generated Stan code parses through rstan when available", {
  skip_if_not_installed("rstan")

  build_components <- getFromNamespace("stan_component_table", "ibis.iSDM")
  build_code <- getFromNamespace("stan_model_code", "ibis.iSDM")

  settings <- Settings$new()
  settings$set("optim_hyperparam", FALSE)

  model <- list(priors = new_waiver())
  model$biodiversity <- list(
    po = list(
      name = "po", type = "poipo", family = "poisson", link = new_waiver(),
      predictors_names = c("p1", "p2"), use_intercept = TRUE
    ),
    pa = list(
      name = "pa", type = "poipa", family = "binomial", link = new_waiver(),
      predictors_names = c("p1", "p2"), use_intercept = TRUE
    )
  )
  model$stan <- list(
    components = build_components(model),
    feature_names = c("p1", "p2"),
    use_dataset_intercepts = TRUE
  )

  po_model <- list(priors = new_waiver())
  po_model$biodiversity <- list(
    po = list(
      name = "po", type = "poipo", family = "poisson", link = new_waiver(),
      predictors_names = c("p1", "p2"), use_intercept = TRUE
    )
  )
  po_model$stan <- list(
    components = build_components(po_model),
    feature_names = c("p1", "p2"),
    use_dataset_intercepts = TRUE
  )

  pa_model <- list(priors = new_waiver())
  pa_model$biodiversity <- list(
    pa = list(
      name = "pa", type = "poipa", family = "binomial", link = new_waiver(),
      predictors_names = c("p1", "p2"), use_intercept = TRUE
    )
  )
  pa_model$stan <- list(
    components = build_components(pa_model),
    feature_names = c("p1", "p2"),
    use_dataset_intercepts = TRUE
  )

  expect_true(rstan::stanc(model_code = build_code(model, settings),
                           model_name = "ibis_stan_multi_test")$status)
  expect_true(rstan::stanc(model_code = build_code(po_model, settings),
                           model_name = "ibis_stan_po_test")$status)
  expect_true(rstan::stanc(model_code = build_code(pa_model, settings),
                           model_name = "ibis_stan_pa_test")$status)
})

test_that("posterior_predict_stanfit orders beta and uses selected component intercept", {
  fake_fit <- structure(
    list(draws = function() {
      posterior::as_draws_array(
        array(
          c(10, 10, 1, 1, 100, 100, 0, 0),
          dim = c(2, 1, 4),
          dimnames = list(NULL, NULL,
                          c("beta[2]", "beta[1]", "Intercept[1]", "Intercept[2]"))
        )
      )
    }),
    class = "CmdStanFit"
  )

  newdata <- data.frame(observed = 0L,
                        p1 = c(1, 0),
                        p2 = c(0, 1))

  out <- posterior_predict_stanfit(
    fake_fit,
    observed ~ 0 + p1 + p2,
    newdata = newdata,
    type = "response",
    link = "logit",
    offset = c(1, 0),
    intercept = "Intercept[2]",
    feature_names = c("p1", "p2")
  )

  expect_true(all(c("mean", "sd", "q05", "q50", "q95", "cv") %in% names(out)))
  expect_equal(out$mean[1], stats::plogis(2), tolerance = 1e-8)
  expect_equal(out$mean[2], stats::plogis(10), tolerance = 1e-8)
  expect_equal(out$sd, c(0, 0), tolerance = 1e-8)
})

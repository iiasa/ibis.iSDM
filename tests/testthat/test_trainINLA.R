context('Train some model with INLA')

# First check that INLA works
test_that('Check that INLA works', {

  skip_if_not_installed('INLA')

  library(INLA)

  # Use test data that comes with INLA
  data(Epil)
  observed <- Epil[1:30, 'y']

  Epil <- rbind(Epil, Epil[1:30, ])
  Epil[1:30, 'y'] <- NA

  # Set up formula and train
  formula = y ~ Trt + Age + V4 + f(Ind, model="iid") + f(rand,model="iid")
  result = inla(formula, family="poisson", data = Epil, control.predictor = list(compute = TRUE, link = 1))

  expect_type(result,'list')
  expect_null(result$waic)
  expect_true(result$ok)
  expect_equal(nrow(result$summary.random$rand), 236)

})


# Train a full distribution model with INLA
test_that('Train a distribution model with INLA', {

  skip_if_not_installed('INLA')

  skip_on_cran()

  # TODO:

})

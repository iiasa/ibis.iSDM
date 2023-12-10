# First check that INLA works
test_that('Create and add priors', {

  # MH: Quick-and-dirty fix for now
  # skip_if_not_installed('INLA')

  # Create list of priors
  p1 <- priors( INLAPrior(variable = 'bias',type = 'normal',
                                      hyper = c(0,1e6) ) )
  # Check whether all information is set and can be queried
  expect_s3_class(p1,'PriorList')
  expect_equal(p1$length(),1)
  expect_equal(as.character(p1$varnames()),'bias')
  expect_equal(as.character(p1$classes()), 'INLAPrior')
  expect_type(p1$exists('bias'),'character')
  expect_length(p1$get('bias'), 2)
  expect_vector(p1$get('bias'), c(0,1e6))
  expect_equal(p1$priors[[1]]$get(), p1$get('bias'))

  # Check empty priors
  expect_error(priors( INLAPrior(variable = '', type = 'normal') ))

  # Linear not supported anymore
  expect_error(INLAPrior(variable = "test",type = "linear"))

  # Now add another prior
  new <- INLAPrior(variable = 'forest', type = 'normal', hyper = c(2, 0.5))

  # Combine priors in list before adding
  p1$add(new)
  expect_length(p1$varnames(),2)
  expect_vector(p1$varnames(),c('bias','forest'))
  expect_vector(p1$get('bias'), c(0,1e6))
  expect_vector(p1$get('forest'), c(2,0.5))
  expect_null(p1$get('dummy'))

  # --- #
  # Get set some priors
  p <- INLAPrior(variable = "CLC3_211_mean_50km",
                 type = "normal",
                 hyper = c(2, 1000) # Precision priors, thus larger sigmas indicate higher precision
  )
  pp <- priors(p)
  expect_equal(pp$get("CLC3_211_mean_50km"), c(2,1000))
  # Emulating some hacky code in inlabru module to check that it works
  vn <- pp$varnames()[which(pp$varnames() == "CLC3_211_mean_50km")]
  ty <- pp$types()[names(vn)]
  expect_vector(vn, "CLC3_211_mean_50km")

  # --- #
  # Instead of adding, combine priors
  new1 <- INLAPrior(variable = 'forest', type = 'normal', hyper = c(2, 0.5))
  new2 <- INLAPrior(variable = 'crops', type = 'normal', hyper = c(-2, 200))

  test1 <- priors(new1,new1) # Duplicated
  expect_length(test1$varnames(), 1)
  expect_equal(test1$length(),1)

  test2 <- priors(new1, new2)
  expect_length(test2$varnames(), 2)
  expect_vector(test2$get('crops'), c(-2,200))
  expect_vector(test2$get('forest'), c(2,0.5))

  # Remove prior
  test2$rm( test2$exists('crops') )
  expect_equal(test2$length(),1)
  expect_null(test2$get('crops'))

  # Set as list
  test3 <- priors(list(new1,new2))
  expect_vector(test3$varnames(),c('forest','crops'))
  expect_vector(test3$ids())
  expect_equal(test3$length(),2)
  expect_vector(test3$get('forest'), c(2,0.5))

  # Or differently set directly
  test3b <- priors(new1,new2)
  expect_vector(test3$varnames(),c('forest','crops'))
  expect_vector(test3$get('crops'), c(-2,200))

  # Combination
  test4 <- priors(new1)
  expect_error(test4$combine(new2))
  test4$combine(priors(new2))
  expect_equal(test4$length(),2)

  expect_vector(test4$get('crops'), c(-2,200))
  expect_vector(test4$get('forest'), c(2,0.5))

  # --- #
  # Combine BREG priors with same name but different values
  new1 <- BREGPrior(variable = 'forest', hyper = 1, ip = 0.5)
  new2 <- BREGPrior(variable = 'forest', hyper = 0, ip = 1)
  pp <- invisible( suppressMessages(priors(new1,new2)) )
  expect_equal(pp$length(), 1)

  # Combine single prior with lists of priors
  new1 <- BREGPrior(variable = 'forest', hyper = 1, ip = 0.5)
  new2 <- BREGPriors(variable = c("shrubs", "cropland"), hyper = 0, ip = 1)
  pp <- invisible( suppressMessages(priors(new1,new2)) )
  expect_equal(pp$length(), 3)

  # Combine two lists of priors
  new1 <- BREGPriors(variable = c('forest', "secondary"), hyper = 1, ip = 0.5)
  new2 <- BREGPriors(variable = c("shrubs", "cropland"), hyper = 0, ip = 1)
  pp <- invisible( suppressMessages(priors(new1,new2)) )
  expect_equal(pp$length(), 4)

  # Add two duplicate but different prior.
  # Default behaviour is to take the last one and raise warning
  new3 <- INLAPrior(variable = 'crops', type = 'normal', hyper = c(2, 0.5))
  new4 <- INLAPrior(variable = 'crops', type = 'normal', hyper = c(-2, 200))
  invisible( suppressMessages( priors(new3, new4) ) )
  suppressMessages( suppressWarnings(pp <- priors(new3, new4)) )
  expect_equal(pp$length(),1)
  expect_equal(pp$get('crops'),c(-2,200))

  # ------- #
  # Testing spde priors and duplicates
  spde1 <- INLAPrior(variable = 'spde', type = 'prior.range', hyper = c(1, 0.5))
  spde2 <- INLAPrior(variable = 'spde', type = 'prior.sigma', hyper = c(1, 0.05))
  suppressMessages( ss <- priors(spde1,spde2) )
  expect_equal(ss$length(),2)
  # Return hyper with defined variable
  expect_vector(ss$get('spde','prior.range'), c(1, 0.5))

  # Two duplicate INLA priors
  new5 <- INLAPrior(variable = 'crops', type = 'normal', hyper = c(2, 0.5))
  new6 <- INLAPrior(variable = 'crops', type = 'normal', hyper = c(-2, 200))
  suppressMessages(nn <- priors(new5,new6))
  expect_vector(nn$get('crops'), c(-2,200))

  # Two gdb priors that are duplicate
  gdb1 <- GDBPrior(variable = 'bias',hyper = 'positive')
  gdb2 <- GDBPrior(variable = 'bias',hyper = 'negative')
  suppressMessages(gg <- priors(gdb1,gdb2))
  expect_equal(gg$length(),1)
  expect_equal(gg$get('bias'),'negative') # Should be last added one

  # ------- #
  # GDB priors just for reference
  expect_error(GDBPrior('bias','oscillating'))
  pg <- priors(GDBPrior('bias','positive'))
  expect_equal(pg$length(),1)
  expect_equal(pg$get('bias'),'positive')

  # --- #
  # Add INLA priors of different types
  pp1 <- INLAPrior(variable = "bias",type = "normal",hyper = c(0,1))
  pp2 <- INLAPrior(variable = "bias",type = "gaussian",hyper = c(0,1))
  suppressMessages(pp <- priors(pp1,pp2))
  expect_equal(pp$length(), 1)
  expect_equal(pp$types() |> as.character(), "gaussian")

  # --- #
  # Different type combinations
  p1 <- INLAPrior(variable = "bias",type = "normal",hyper = c(0,1))
  p2 <- INLAPrior(variable = "forest",type = "clinear",hyper = c(0,Inf))
  pp1 <- priors(p1,p2)
  p3 <- INLAPrior(variable = "bias",type = "clinear",hyper = c(Inf,0))
  p4 <- INLAPrior(variable = "forest",type = "normal",hyper = c(0,100))
  pp2 <- priors(p3,p4)
  # Combine priors
  expect_invisible(pp1$combine(pp2))
  expect_equal(pp1$get("forest"), c(0,100))

  # --- #
  # Check that name sanitation works.
  pg <- priors(GDBPrior('FunckyName--Test__f','positive'))
  expect_equal(as.character(pg$varnames()), "FunckyName_Test_f")


})

test_that('Add and modify priors to existing object', {

  skip_if_not_installed('INLA')
  skip_on_cran()
  skip_on_ci()

  suppressWarnings( requireNamespace("terra", quietly = TRUE) )
  options("ibis.setupmessages" = FALSE)

  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))
  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE), 'points',quiet = TRUE)
  ll <- list.files(system.file('extdata/predictors/',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)

  predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Create list of priors
  pp <- priors( INLAPrior(variable = 'CLC3_132_mean_50km',type = 'normal',
                          hyper = c(2,1e6) ) )

  # Define a model
  invisible(
    x <- distribution(background) |>
    add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed', name = 'Virtual points') |>
    add_predictors(predictors[[c('slope_mean_50km','bio01_mean_50km','CLC3_132_mean_50km')]],
                   transform = 'none',derivates = 'none') |>
    engine_inla(
      max.edge = c(.5, 3),
      offset = c(0.5, 1),
      cutoff = 0.5,
      proj_stepsize = 1
    )
  )
  expect_s3_class(x$priors,'Waiver')
  expect_s3_class(x$get_priors(),'Waiver')

  # Add priors to it
  x <- x |>  add_priors(priors = pp)

  expect_s3_class( x$get_priors(), 'PriorList')
  expect_vector(x$get_prior_variables(), "CLC3_132_mean_50km" )

  # Remove priors from it
  invisible( x |> rm_priors() )
  expect_s3_class( x$get_priors(), 'PriorList')
  x <- x |> rm_priors()
  expect_s3_class(x$priors,'Waiver')
  expect_s3_class(x$get_priors(),'Waiver')

  # Add duplicated priors to it with the same name
  p1 <- INLAPrior(variable = 'CLC3_132_mean_50km',type = 'normal',
                          hyper = c(2,1e6) )
  p2 <- INLAPrior(variable = 'CLC3_132_mean_50km',type = 'clinear',
                  hyper = c(0,Inf) )

  suppressMessages( expect_equal( priors(p1,p2)$length(), 1 ) )
  pp <- priors(p1,p2)
  expect_vector(pp$get('CLC3_132_mean_50km'), c(0, Inf))

  # With BREG engine
  suppressWarnings( suppressMessages(x <- x |> engine_breg()) )
  p1 <- BREGPrior("test", hyper = 0.5, ip = 1)
  p2 <- BREGPrior("test", hyper = 5, ip = .2)
  expect_equal(priors(p1,p2)$length(),1)

  x <- x |> add_priors(priors(p1,p2))
  expect_equal(x$get_prior_variables() |> as.character(), "test")
  x <- x |> add_priors(priors(
    BREGPrior("test2", hyper = 5, ip = .2)
  ))
  expect_true("test2" %in% x$get_prior_variables() )

  # Combine priors with the same type
  x <- x |> rm_priors()
  expect_s3_class(x, "BiodiversityDistribution")
  p1 <- INLAPrior("test", type = "normal", hyper = c(0, 0.05))
  p2 <- INLAPrior("test", type = "clinear", hyper = c(0, Inf))
  x <- x |> add_priors( priors(p1,p2) )
  # Reverse and add again
  x <- x |> add_priors( priors(p2, p1) )
  expect_equal(x$priors$length(), 1)

  # Check SPDE priors in multiple versions
  x <- x |> rm_priors()
  spde1 <- INLAPrior(variable = 'spde', type = 'prior.range', hyper = c(1, 0.5))
  spde2 <- INLAPrior(variable = 'spde', type = 'prior.sigma', hyper = c(1, 0.05))
  ss <- priors(spde1, spde2)
  x <- x |> add_priors(ss)
  expect_equal(x$priors$length(), 2)
  # Now reverse and add again. Should still be 2 priors
  x <- x |> add_priors(priors(spde2, spde1))
  expect_equal(x$priors$length(), 2)

  # Summarize priors
  pp <- x$get_priors()
  expect_s3_class(summary(pp), "data.frame")
})

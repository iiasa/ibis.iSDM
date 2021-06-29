context('Prior settings')

# First check that INLA works
test_that('Create and add priors', {

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

  # Now add another prior
  new <- INLAPrior(variable = 'forest', type = 'normal', hyper = c(2, 0.5))
  p1$add(new)

  expect_length(p1$varnames(),2)
  expect_vector(p1$varnames(),c('bias','forest'))
  expect_vector(p1$get('bias'), c(0,1e6))
  expect_vector(p1$get('forest'), c(2,0.5))
  expect_null(p1$get('dummy'))

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

  # Combination
  test3 <- priors(new1)
  expect_error(test3$combine(new2))
  test3$combine(priors(new2))
  expect_equal(test3$length(),2)

  expect_vector(test3$get('crops'), c(-2,200))
  expect_vector(test3$get('forest'), c(2,0.5))

  # ------- #
  # GDB priors just for reference
  expect_error(GDBPrior('bias','oscillating'))
  pg <- priors(GDBPrior('bias','positive'))
  expect_equal(pg$length(),1)
  expect_equal(pg$get('bias'),'positive')
})

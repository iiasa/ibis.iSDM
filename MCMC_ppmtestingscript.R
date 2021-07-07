# Libraries
library(spatstat)
library(sf)
library(sp)
library(maptools)
library(raster)
library(fields)
library(rstan)
library(tidyverse)
library(RandomFields)
library(bayesplot)

# Simulate realization of a cox process
genDat_cox <- function(b0, b1, dim, noise_mean = NULL, noise_sd = NULL, plotdat = TRUE){

  # Define the window of interest
  win <- owin(c(0,dim[1]), c(0,dim[2]))

  # set number of pixels to simulate an environmental covariate
  spatstat.options(npixel=c(dim[1],dim[2]))

  y0 <- seq(win$yrange[1], win$yrange[2],
            length=spatstat.options()$npixel[2])
  x0 <- seq(win$xrange[1], win$xrange[2],
            length=spatstat.options()$npixel[1])
  multiplier <- 1/dim[2]

  # Make the environmental covariate
  gridcov <- outer(x0,y0, function (x,y) multiplier*y + 0*x)

  # Set the parameter values
  beta0 <- b0
  beta1 <- b1

  if(!is.null(noise_mean) && !is.null(noise_sd)){
    noise_mean <- noise_mean
    noise_sd <- noise_sd
  }

  else{
    noise_mean = 0
    noise_sd = 1
  }
  # Create 'im' objects for simulating the point process
  # First we create a random field (just noise), then the intensity
  # field made of our linear predictors and we sum up the two images
  # to get the intensity of the point process
  noise <- rnoise(rnorm, mean = noise_mean, sd = noise_sd, w = win)
  linear <- im(b0 + b1*gridcov, xrange = c(0,20), yrange = c(0,20))
  intensity <- noise + linear

  # Simulate the point pattern
  pp <- rpoispp(exp(intensity), xcol=x0, yrow=y0)
  qcounts <- quadratcount(pp, ny=dim[1], nx=dim[2])
  dens <- density(pp)
  Lambda <- as.vector(t(qcounts))

  if(plotdat == TRUE){
    par(mfrow=c(2,2), mar=c(2,2,1,1), mgp=c(1,0.5,0))
    plot(noise, main = 'White noise')
    plot(im(gridcov), main = 'Covariate')
    plot(intensity, main = 'log Intensity')
    plot(dens, main = 'Intensity of the point pattern')
  }
  # Return a list with which I can play with
  return(list(Lambda = Lambda, pp = pp, gridcov = gridcov))
}

# Set a seed
set.seed(123)

# We now have a double stochastic process where the intensity is random
b0 <- 2
b1 <- 3
dim <- c(20,20)
noise_mean <- 1
noise_sd <- 0.5

# Generate data
pp <- genDat_cox(b0, b1, dim, noise_mean, noise_sd,plotdat = TRUE)


# ------ #
# Cox stan model -----
# Recover those parameters in stan!

cox_stan <- '
// Fit a Cox process in Stan

data{
  int<lower = 1> n;
  vector[n] x;
  int<lower = 0> y[n];
}
parameters{
  real beta0;
  real beta1;
  real<lower = 0, upper = 5> sigma_noise;
  vector[n] noise;
}
transformed parameters{
}
model{
  //priors
  target += normal_lpdf(beta0 | 0,5);
  target += normal_lpdf(beta1 | 0,10);
  target += uniform_lpdf(sigma_noise | 0,1);

  // Prior for the noise
  target += normal_lpdf(noise | 0, sigma_noise);

  // likelihood
  target += poisson_log_lpmf(y | beta0 + beta1 * x + noise);
}
generated quantities{
  vector[n] lambda_rep;
  lambda_rep = exp(beta0 + beta1 * x + noise);
}
'
# Prepare the data list for Stan
stan_data = list(n = length(pp$Lambda), x = as.vector(t(pp$gridcov)), y = pp$Lambda)

# Fit the model
fit_stan_cox <- stan(model_code = cox_stan, data = stan_data,
                     warmup = 500, iter = 2000, chains = 3)

# Take a look at the model output
print(fit_stan_cox, pars = c('beta0', 'beta1', 'sigma_noise'))

# Get the posterior distribution of the parameters
posterior <- as.array(fit_stan_cox)

# Plot!
mcmc_intervals(posterior,
               pars = c('beta0', 'beta1', 'sigma_noise'),
               prob = 1)

# Retrieve the predicted lambdas from the posterior and calculate their means and standard deviation
lambda_rep <- as.data.frame(rstan::extract(fit_stan_cox)['lambda_rep'])
mean_lambda_rep <- apply(lambda_rep, 2, 'mean')
sd_lambda_rep <- apply(lambda_rep, 2, 'sd')

# Transform the pp object into an sf object so we can count the number of points in each grid cell
pointp <- pp$pp
pp_sp <- as.SpatialPoints.ppp(pointp)
pp_sf <- st_as_sf(pp_sp)

# Create a grid and place the predictions in it
grid <- st_make_grid(pp_sf, n = dim, what = 'polygons') %>%
  st_as_sf() %>%
  mutate(pred = mean_lambda_rep,
         sd = sd_lambda_rep)

# COunt the number of points in each grid cell
grid$real <- lengths(st_intersects(grid, pp_sf))

# Plot the grid
plot(grid["pred"])
plot(grid["sd"])

# LGCP with random field ----

# Simulation of the random field
noise <- rnoise(rnorm, mean = noise_mean, sd = noise_sd)
plot(noise)

noise_smooth <- Smooth(noise, sigma=2, normalise=TRUE, bleed=FALSE)
plot(noise_smooth)

# Simulate realization of a Log-Gaussian Cox process
genDat_lgcp <- function(b0, b1, dim, var, scale, plotdat = TRUE){

  # Define the window of interest
  win <- owin(c(0,dim[1]), c(0,dim[2]))

  # set number of pixels to simulate an environmental covariate
  spatstat.options(npixel=c(dim[1],dim[2]))

  y0 <- seq(win$yrange[1], win$yrange[2],
            length=spatstat.options()$npixel[2])
  x0 <- seq(win$xrange[1], win$xrange[2],
            length=spatstat.options()$npixel[1])
  multiplier <- 1/dim[2]

  # Make the environmental covariate
  gridcov <- outer(x0,y0, function (x,y) multiplier*y + 0*x)

  # Set the parameter values
  beta0 <- b0
  beta1 <- b1
  var <- var
  scale <- scale

  # Simulate the LGCP, here we define the covariance structure as being exponential
  GP <- rLGCP(model="exp",
              mu=im(beta0 + beta1*gridcov, xcol=x0, yrow=y0),
              var=var, scale=scale, win = win)

  # Get the realisation of the LGCP as an sf object - easier to handle
  g <- as.ppp(GP)
  GP_sp <- as.SpatialPoints.ppp(g)
  GP_sf <- st_as_sf(GP_sp)

  # Get the result in a grid
  grid <- st_make_grid(GP_sf, n = dim, what = 'polygons') %>%
    st_as_sf() %>%
    mutate(Lambda = lengths(st_intersects(., GP_sf)),
           cov = as.vector(t(gridcov)))

  if(plotdat == TRUE){
    par(mfrow=c(1,2), mar=c(2,2,1,1), mgp=c(1,0.5,0))
    plot(grid["Lambda"], main = 'Intensity of the point pattern')
  }
  # Return a list with which I can play with
  return(grid)
}

beta0 <- 2
beta1 <- 3
var <- 0.5
scale <- 0.4
dim = c(10,10)

data_lgcp <- genDat_lgcp(beta0, beta1, dim, var, scale)


fit_lgcp0 <- '
// Fit an accurate LGCP in Stan with the Exponential covariance structure
functions{

    matrix GP(matrix x, real sigma_sq, real scale, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sigma_sq + delta;
          for (j in (i + 1):N) {
            K[i, j] = sigma_sq * exp(- x[i,j] / scale );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sigma_sq + delta;
        return K;
    }
}
data{
  int<lower = 1> N;
  vector[N] x;
  int<lower = 0> y[N];
  matrix[N, N] DMat; // Distance matrix
}
parameters{
  real beta0;
  real beta1;
  vector[N] k;

  // GP standard deviation parameters
  real<lower=0> sigma_sq;
  // GP length-scale parameters
  real<lower=0> scale;
}
model{
  matrix[N,N] SIGMA;
  vector[N] mu;

  SIGMA = GP(DMat, sigma_sq, scale, 0.01);
  k ~ multi_normal(rep_vector(0,N), SIGMA);

  //priors for the coefficients
  target += normal_lpdf(beta0 | 0,5);
  target += normal_lpdf(beta1 | 0,10);

  // Prior for the noise
  target += cauchy_lpdf(sigma_sq | 0, 1);
  target += inv_gamma_lpdf(scale | 3.884416, 0.77454);

  // likelihood
    for(i in 1:N){
    mu[i] = beta0 + beta1 * x[i] + k[i];
  }

  target += poisson_log_lpmf(y | mu);
}
generated quantities{
}
'

# To calculate in stan first calculate Dmat:
DMat <- st_distance(st_centroid(data_lgcp), by_element = FALSE)

# Make stan data
stan_data <- list(N = nrow(data_lgcp),
                  x = data_lgcp$cov,
                  y = data_lgcp$Lambda,
                  DMat = DMat)

# Compute the distance matrix
stan_fit0 <- stan(model_code = fit_lgcp0,
                  data = stan_data,
                  chains = 1, warmup = 1000, iter = 5000,
                  control = list(adapt_delta = 0.999, max_treedepth=13))


# Get the posterior of the parameters
draws <- rstan::extract(stan_fit0, pars = c('beta0', 'beta1', 'sigma_sq', 'scale'))

# We make a sequence of distance
dist_seq <- seq(from = min(DMat), to = max(DMat), length.out = 100)

# Compute the mean and the standard deviation of the posterior correlation
post_cov <- sapply(dist_seq,function(x)draws$sigma_sq*exp(-draws$scale*x^2))
post_cov_mu <-apply(post_cov,2,mean)
post_cov_sd <-apply(post_cov,2,sd)

# Make a dataframe and plot
post_df <- tibble(dist = dist_seq,
                  mu = post_cov_mu,
                  sd = post_cov_sd)

ggplot(post_df, aes(x = dist)) +
  geom_line(aes(y = mu), color = "#CD295A", size = 1) +
  geom_ribbon(aes(ymin = mu - sd, ymax = mu + sd), fill = "#38ADAE", alpha = .3) +
  theme_classic() +
  ylab("Covariance") +
  xlab("Distance")

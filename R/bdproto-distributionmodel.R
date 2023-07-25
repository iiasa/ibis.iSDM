#' @include utils.R bdproto.R identifier.R
NULL

#' @export
if (!methods::isClass("DistributionModel")) methods::setOldClass("DistributionModel")
NULL

#' Prototype for the trained Model object
#'
#' All trained Models should inherit the options here
#'
#' @name DistributionModel-class
#' @aliases DistributionModel
#' @family bdproto
#' @keywords bdproto
NULL

#' @export
DistributionModel <- bdproto(
  "DistributionModel",
  id = character(), # An id for any trained model
  model = list(),
  fits = list(), # List of fits with data
  # Print message with summary of model
  print = function(self) {
    # TODO: Have a lot more information in here and to be prettified

    # Check whether prediction exists and number of layers
    has_prediction <- "prediction" %in% self$show_rasters()
    # Check whether threshold has been calculated
    has_threshold <- grep('threshold',self$show_rasters(),value = TRUE)[1]

    # FIXME: Have engine-specific code moved to engine
    if(inherits(self, 'INLA-Model') || inherits(self, 'INLABRU-Model') ){
      if( length( self$fits ) != 0 ){
        # Get strongest effects
        ms <- subset(tidy_inla_summary(self$get_data('fit_best')),
                     select = c('variable', 'mean'))
        ms <- ms[order(ms$mean,decreasing = TRUE),] # Sort

        message(paste0(
          'Trained ',class(self)[1],' (',self$show(),')',
          '\n  \033[2mStrongest summary effects:\033[22m',
          '\n     \033[34mPositive:\033[39m ', name_atomic(ms$variable[ms$mean>0]),
          '\n     \033[31mNegative:\033[39m ', name_atomic(ms$variable[ms$mean<0]),
          ifelse(has_prediction,
                 paste0("\n  Prediction fitted: ",text_green("yes")),
                 ""),
          ifelse(!is.na(has_threshold),
                 paste0("\n  Threshold created: ",text_green("yes")),
                 "")
        ))
      }
    } else if( inherits(self, 'GDB-Model') ) {

        # Get Variable importance
        vi <- mboost::varimp(
          self$get_data('fit_best')
        )
        vi <- sort( vi[which(vi>0)],decreasing = TRUE )

        message(paste0(
          'Trained ',class(self)[1],' (',self$show(),')',
          '\n  \033[2mStrongest effects:\033[22m',
          '\n     ', name_atomic(names(vi)),
          ifelse(has_prediction,
                 paste0("\n  Prediction fitted: ",text_green("yes")),
                 ""),
          ifelse(!is.na(has_threshold),
                 paste0("\n  Threshold created: ",text_green("yes")),
                 "")
        ))
    } else if( inherits(self, 'BART-Model') ) {
      # Calculate variable importance from the posterior trees
      vi <- varimp.bart(self$get_data('fit_best'))

      message(paste0(
        'Trained ',class(self)[1],' (',self$show(),')',
        '\n  \033[2mStrongest effects:\033[22m',
        '\n     ', name_atomic(vi$names),
        ifelse(has_prediction,
               paste0("\n  Prediction fitted: ",text_green("yes")),
               ""),
        ifelse(!is.na(has_threshold),
               paste0("\n  Threshold created: ",text_green("yes")),
               "")
      ))
    } else if( inherits(self, 'STAN-Model') ) {
      # Calculate variable importance from the posterior
      vi <- rstan::summary(self$get_data('fit_best'))$summary |> as.data.frame() |>
        tibble::rownames_to_column(var = "parameter") |> as.data.frame()
      # Get beta coefficients only
      vi <- vi[grep("beta", vi$parameter,ignore.case = TRUE),]

      # Get variable names from model object
      # FIXME: This might not work for all possible modelling objects. For instance
      model <- self$model
      assertthat::assert_that(nrow(vi) == length(model$predictors_names),
                              length(vi$parameter) == length(model$predictors_names))
      vi$parameter <- model$predictors_names

      vi <- vi[order(abs(vi$mean),decreasing = TRUE),]
      message(paste0(
        'Trained ',class(self)[1],' (',self$show(),')',
        '\n  \033[2mStrongest summary effects:\033[22m',
        '\n     \033[34mPositive:\033[39m ', name_atomic(vi$parameter[vi$mean>0]),
        '\n     \033[31mNegative:\033[39m ', name_atomic(vi$parameter[vi$mean<0])
      ))
    } else if( inherits(self, 'XGBOOST-Model') ) {
      vi <- xgboost::xgb.importance(model = self$get_data('fit_best'),)

      message(paste0(
        'Trained ',class(self)[1],' (',self$show(),')',
        '\n  \033[2mStrongest effects:\033[22m',
        '\n     ', name_atomic(vi$Feature),
        ifelse(has_prediction,
               paste0("\n  Prediction fitted: ",text_green("yes")),
               ""),
        ifelse(!is.na(has_threshold),
               paste0("\n  Threshold created: ",text_green("yes")),
               "")
      ))
    } else if( inherits(self, 'BREG-Model') ) {
      obj <- self$get_data('fit_best')
      # Summarize the beta coefficients from the posterior
      ms <- posterior::summarise_draws(obj$beta) |>
        subset(select = c('variable', 'mean'))
      # Reorder
      ms <- ms[order(ms$mean,decreasing = TRUE),] # Sort

      message(paste0(
        'Trained ',class(self)[1],' (',self$show(),')',
        '\n  \033[2mStrongest summary effects:\033[22m',
        '\n     \033[34mPositive:\033[39m ', name_atomic(ms$variable[ms$mean>0]),
        '\n     \033[31mNegative:\033[39m ', name_atomic(ms$variable[ms$mean<0]),
        ifelse(has_prediction,
               paste0("\n  Prediction fitted: ",text_green("yes")),
               ""),
        ifelse(!is.na(has_threshold),
               paste0("\n  Threshold created: ",text_green("yes")),
               "")
      ))
    } else if(inherits(self, 'GLMNET-Model')) {
      obj <- self$get_data('fit_best')

      # Summarise coefficients within 1 standard deviation
      ms <- tidy_glmnet_summary(obj)

      message(paste0(
        'Trained ',class(self)[1],' (',self$show(),')',
        '\n  \033[1mStrongest summary effects:\033[22m',
        '\n     \033[34mPositive:\033[39m ', name_atomic(ms$variable[ms$mean>0]),
        '\n     \033[31mNegative:\033[39m ', name_atomic(ms$variable[ms$mean<0]),
        ifelse(has_prediction,
               paste0("\n  Prediction fitted: ",text_green("yes")),
               ""),
        ifelse(!is.na(has_threshold),
               paste0("\n  Threshold created: ",text_green("yes")),
               "")
      ))

    } else {
      message(paste0(
        'Trained distribution model (',self$show(),')',
        text_red('\n     No fitted model found!'),
        ifelse(has_prediction,
               paste0("\n  Prediction fitted: ",text_green("yes")),
               ""),
        ifelse(!is.na(has_threshold),
               paste0("\n  Threshold created: ",text_green("yes")),
               "")
      ))
    }
  },
  # Show the name of the Model
  show = function(self) {
    self$model$runname
  },
  # Plot the prediction
  plot = function(self, what = 'mean'){
    if( length( self$fits ) != 0 && !is.null( self$fits$prediction ) ){
      pred <- self$get_data('prediction')
      assertthat::assert_that(is.Raster(pred))
      # Check if median is requested but not present, change to q50
      if(what == "median" && !(what %in% names(pred))) { what <- "q50" }

      # Match argument
      what <- match.arg(what, names(pred), several.ok = FALSE)
      assertthat::assert_that( what %in% names(pred),msg = paste0('Prediction type not found. Available: ', paste0(names(pred),collapse = '|')))
      terra::plot(pred[[what]],
           main = paste0(self$model$runname, ' prediction (',what,')'),
           box = FALSE,
           axes = TRUE,
           colNA = NA, col = ibis_colours[['sdm_colour']]
           )
    } else {
      message(
        paste0('No model predictions found.')
      )
    }
  },
  # Plot threshold
  plot_threshold = function(self, what = 1){
    assertthat::assert_that(is.numeric(what) || is.character(what))
    # Determines whether a threshold exists and plots it
    rl <- self$show_rasters()
    if(length(grep('threshold',rl))>0){

      # Get stack of computed thresholds
      ras <- self$get_data( grep('threshold', rl, value = TRUE) )[[what]]
      suppressWarnings(
        ras <- terra::droplevels(ras)
      )
      # Get colour palette
      format <- attr(ras[[1]], 'format') # Format attribute
      if(is.null(format)) format = "binary"
      if(format == "normalize"){
        col <- colorRampPalette(c("grey","#EB072F","#FFE900","#5A94DD","black"))(100)
      } else if(format == "percentile") {
        col <- colorRampPalette(c("grey","#EB072F","#FFE900","#5A94DD","black"))(length(unique(ras)[,1]))
      } else {
        # Binary
        col <- c("grey", "black")
      }
      terra::plot(ras,
                   box = FALSE,
                   axes = TRUE,
                   colNA = NA, col = col
      )
    } else {
      message("No computed threshold was found!")
      invisible()
    }
  },
  # Show model run time if settings exist
  show_duration = function(self){
    if(!is.Waiver(self$settings)) self$settings$duration()
  },
  # Get effects or importance tables from model
  summary = function(self, obj = 'fit_best'){
    # Distinguishing between model types
    if(inherits(self, 'GDB-Model')){
      clean_mboost_summary( self$get_data(obj) )
    } else if(inherits(self, 'INLA-Model') || inherits(self, 'INLABRU-Model')){
      tidy_inla_summary(self$get_data(obj))
    } else if(inherits(self, 'BART-Model')){
      # Number of times each variable is used by a tree split
      # Tends to become less informative with higher numbers of splits
      varimp.bart(self$get_data(obj)) |> tibble::remove_rownames()
    } else if(inherits(self, 'STAN-Model')){
      vi <- rstan::summary(self$get_data(obj))$summary |> as.data.frame() |>
        tibble::rownames_to_column(var = "parameter") |> as.data.frame()
      # Get beta coefficients only
      vi <- vi[grep("beta", vi$parameter,ignore.case = TRUE),]
      # FIXME: This might not work for all possible modelling objects. For instance
      model <- self$model
      assertthat::assert_that(nrow(vi) == length(model$predictors_names),
                              length(vi$parameter) == length(model$predictors_names))
      vi$parameter <- model$predictors_names
      names(vi) <- make.names(names(vi))
      return( tibble::as_tibble( vi ) )
    } else if(inherits(self, 'BREG-Model')){
      posterior::summarise_draws(self$get_data(obj)$beta)
    } else if(inherits(self, "XGBOOST-Model")){
      xgboost::xgb.importance(model = self$get_data(obj))
    } else if(inherits(self, 'GLMNET-Model')){
      tidy_glmnet_summary(self$get_data(obj))
    }
  },
  # Model convergence check
  has_converged = function(self){
    new_waiver()
  },
  # Dummy residual function
  get_residuals = function(self){
    new_waiver()
  },
  # Dummy partial response calculation. To be overwritten per engine
  partial = function(self){
    new_waiver()
  },
  # Dummy spartial response calculation. To be overwritten per engine
  spartial = function(self){
    new_waiver()
  },
  # Generic plotting function for effect plots
  effects = function(self, x = 'fit_best', what = 'fixed', ...){
    assertthat::assert_that(is.character(what))
    if(inherits(self, 'GDB-Model')){
      # How many effects
      n <- length( stats::coef( self$get_data(x) ))
      # Use the base plotting
      par.ori <- graphics::par(no.readonly = TRUE)
      graphics::par(mfrow = c(ceiling(n/3),3))

      mboost:::plot.mboost(x = self$get_data(x),
                           type = 'b',cex.axis=1.5, cex.lab=1.5)

      graphics::par(par.ori)#dev.off()
    } else if(inherits(self, 'INLA-Model')) {
      plot_inla_marginals(self$get_data(x),what = what)
    } else if(inherits(self, 'GLMNET-Model')) {
      if(what == "fixed"){
        glmnet:::plot.glmnet(self$get_data(x)$glmnet.fit, xvar = "lambda") # Deviance explained
      } else{ plot(self$get_data(x)) }
    } else if(inherits(self, 'STAN-Model')) {
      # Get true beta parameters
      ra <- grep("beta", names(self$get_data(x)),value = TRUE) # Get range
      rstan::stan_plot(self$get_data(x), pars = ra)
    } else if(inherits(self, 'INLABRU-Model')) {
      # Use inlabru effect plot
      ggplot2::ggplot() +
        inlabru::gg(self$get_data(x)$summary.fixed, bar = TRUE)
    } else if(inherits(self, 'BART-Model')){
      message('Calculating partial dependence plots')
      self$partial(self$get_data(x), x.vars = what, ...)
    } else if(inherits(self, 'BREG-Model')){
      obj <- self$get_data(x)
      if(what == "fixed") what <- "coefficients"
      what <- match.arg(what, choices = c("coefficients", "scaled.coefficients","residuals",
                                           "size", "fit", "help", "inclusion"), several.ok = FALSE)
      if( length( grep("poisson", obj$call) ) > 0 ){
        BoomSpikeSlab::plot.poisson.spike(obj, y = what)
      } else if( length( grep("binomial", obj$call) ) > 0 ){
        BoomSpikeSlab::plot.logit.spike(obj, y = what)
      } else {
        BoomSpikeSlab::plot.lm.spike(obj, y = what)
      }
    } else if(inherits(self, "XGBOOST-Model")){
      # Check whether linear model was fitted, otherwise plot tree
      if( self$settings$get("only_linear") ){
        vi <- self$summary(x)
        xgboost::xgb.ggplot.importance(vi)
      } else {
        obj <- self$get_data(x)
        xgboost::xgb.plot.multi.trees(obj)
      }
    } else {
      self$partial(self$get_data(x), x.vars = NULL)
    }
  },
  # Get equation
  get_equation = function(self){
    self$get_data("fit_best_equation")
  },
  # Get specific fit from this Model
  get_data = function(self, x = "prediction") {
    if (!x %in% names(self$fits))
      return(new_waiver())
    return(self$fits[[x]])
  },
  # Set fit for this Model
  set_data = function(self, x, value) {
    # Get biodiversity dataset collection
    ff <- self$fits
    # Set the object
    ff[[x]] <- value
    bdproto(NULL, self, fits = ff )
  },
  # Get the threshold value if calculated
  get_thresholdvalue = function(self){
    # Determines whether a threshold exists and plots it
    rl <- self$show_rasters()
    if(length(grep('threshold',rl))==0) return( new_waiver() )

    # Get the thresholded layer and return the respective attribute
    obj <- self$get_data( grep('threshold',rl,value = TRUE) )
    assertthat::assert_that(assertthat::has_attr(obj, "threshold"))
    return(
      attr(obj, "threshold")
    )
  },
  # List all rasters in object
  show_rasters = function(self){
    rn <- names(self$fits)
    rn <- rn[ which( sapply(rn, function(x) is.Raster(self$get_data(x)) ) ) ]
    return(rn)
  },
  # Get projection
  get_projection = function(self){
    sf::st_crs(self$model$background)
  },
  # Get resolution
  get_resolution = function(self){
    if(!is.Waiver(self$get_data())){
      terra::res( self$get_data() )
    } else {
      # Try to get it from the modelling object
      self$model$predictors_object$get_resolution()
    }
  },
  # Remove calculated thresholds
  rm_threshold = function(self){
    rl <- self$show_rasters()
    if(length(grep('threshold',rl))>0){
      for(val in grep('threshold',rl,value = TRUE)){
        self$fits[[val]] <- NULL
      }
    }
    invisible()
  },
  # Calculate a suitability index
  calc_suitabilityindex = function(self, method = "normalize"){
    assertthat::assert_that(
      is.character(method),
      is.Raster(self$get_data())
    )
    method <- match.arg(method, c("normalize", "reltotal"), several.ok = FALSE)

    # Get the raster of the mean prediction
    ras <- self$get_data()[["mean"]]
    if(method == "normalize"){
      out <- predictor_transform(ras, option = "norm")
    } else {
      out <- ras / terra::global(ras,"sum", na.rm = TRUE)[,1]
    }
    return(out)
  },
  # Save object
  save = function(self, fname, type = 'gtif', dt = 'FLT4S'){
    assertthat::assert_that(
      is.character(fname),
      type %in% c('gtif','gtiff','tif','nc','ncdf'),
      'fits' %in% self$ls(),
      dt %in% c('LOG1S','INT1S','INT1U','INT2S','INT2U','INT4S','INT4U','FLT4S','FLT8S')
    )
    type <- tolower(type)

    # Get raster file in fitted object
    cl <- sapply(self$fits, class)
    ras <- self$fits[[grep('raster', cl,ignore.case = T)]]

    # Check that no-data value is not present in ras
    assertthat::assert_that(any(!terra::global(ras, "min", na.rm = TRUE)[,1] <= -9999),
                            msg = 'No data value -9999 is potentially in prediction!')

    if(file.exists(fname)) warning('Overwritting existing file...')
    if(type %in% c('gtif','gtiff','tif')){
      # Save as geotiff
      writeGeoTiff(ras, fname = fname, dt = dt)
    } else if(type %in% c('nc','ncdf')) {
      # Save as netcdf
      # TODO: Potentially change the unit descriptions
      writeNetCDF(ras, fname = fname, varName = 'iSDM prediction', varUnit = "",varLong = "")
    }
    invisible()
  }
)

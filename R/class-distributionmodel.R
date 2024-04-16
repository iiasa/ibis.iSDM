if (!methods::isClass("DistributionModel")) methods::setOldClass("DistributionModel")

#' Class for the trained Model object
#'
#' @description
#' All trained Models inherit the options here plus any additional ones defined
#' by the engine and inference.
#'
#' @keywords classes
#'
#' @name DistributionModel-class
NULL

#' @rdname DistributionModel-class
#' @export
DistributionModel <- R6::R6Class(
  "DistributionModel",
  public = list(
    #' @field id A character id for any trained model
    #' @field name A description of the model as [`character`].
    #' @field model A [`list`] containing all input datasets and parameters to the model.
    #' @field settings A [`Settings-class`] object with information on inference.
    #' @field fits A [`list`] containing the prediction and fitted model.
    #' @field .internals A [`list`] containing previous fitted models.
    id = character(),
    name = character(),
    model = list(),
    settings = new_waiver(),
    fits = list(),
    .internals = new_waiver(),

    #' @description
    #' Initializes the object and creates an empty list
    #' @param name A description of the model as [`character`].
    #' @return NULL
    initialize = function(name){
      assertthat::assert_that(
        is.character(name)
      )
      self$name <- name
    },

    #' @description
    #' Return the name of the model
    #' @return A [`character`] with the model name used.
    get_name = function(){
      return( self$name )
    },

    #' @description
    #' Print the names and summarizes the model within
    #' @note
    #' Could be further pretified and commands outsourced.
    #' @return A message on screen
    print = function() {
      # Check whether prediction exists and number of layers
      has_prediction <- "prediction" %in% self$show_rasters()
      # Check whether threshold has been calculated
      has_threshold <- grep('threshold',self$show_rasters(),value = TRUE)[1]

      # Get model
      obj <- self$get_data('fit_best')

      # FIXME: Have engine-specific code moved to engine
      if( self$get_name() == 'INLA-Model' || self$get_name() == 'INLABRU-Model'){
        if( length( self$fits ) != 0 ){
          # Get strongest effects
          ms <- subset(tidy_inla_summary(obj),
                       select = c('variable', 'mean'))
          ms <- ms[order(ms$mean,decreasing = TRUE),] # Sort

          message(paste0(
            'Trained ',self$name,' (',self$show(),')',
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
      } else if( self$get_name() == 'GDB-Model' ) {

        # Get Variable importance
        vi <- mboost::varimp(obj)
        vi <- sort( vi[which(vi>0)],decreasing = TRUE )

        message(paste0(
          'Trained ',self$name,' (',self$show(),')',
          '\n  \033[2mStrongest effects:\033[22m',
          '\n     ', name_atomic(names(vi)),
          ifelse(has_prediction,
                 paste0("\n  Prediction fitted: ",text_green("yes")),
                 ""),
          ifelse(!is.na(has_threshold),
                 paste0("\n  Threshold created: ",text_green("yes")),
                 "")
        ))
      } else if( self$get_name() == 'BART-Model' ) {
        # Calculate variable importance from the posterior trees
        vi <- varimp.bart(obj)

        message(paste0(
          'Trained ',self$name,' (',self$show(),')',
          '\n  \033[2mStrongest effects:\033[22m',
          '\n     ', name_atomic(vi$names),
          ifelse(has_prediction,
                 paste0("\n  Prediction fitted: ",text_green("yes")),
                 ""),
          ifelse(!is.na(has_threshold),
                 paste0("\n  Threshold created: ",text_green("yes")),
                 "")
        ))
      } else if( self$get_name() == 'STAN-Model' ) {
        # Calculate variable importance from the posterior
        vi <- rstan::summary(obj)$summary |> as.data.frame() |>
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
          'Trained ',self$name,' (',self$show(),')',
          '\n  \033[2mStrongest summary effects:\033[22m',
          '\n     \033[34mPositive:\033[39m ', name_atomic(vi$parameter[vi$mean>0]),
          '\n     \033[31mNegative:\033[39m ', name_atomic(vi$parameter[vi$mean<0])
        ))
      } else if( self$get_name() == 'XGBOOST-Model' ) {
        vi <- xgboost::xgb.importance(model = obj)

        message(paste0(
          'Trained ',self$name,' (',self$show(),')',
          '\n  \033[2mStrongest effects:\033[22m',
          '\n     ', name_atomic(vi$Feature),
          ifelse(has_prediction,
                 paste0("\n  Prediction fitted: ",text_green("yes")),
                 ""),
          ifelse(!is.na(has_threshold),
                 paste0("\n  Threshold created: ",text_green("yes")),
                 "")
        ))
      } else if( self$get_name() == 'BREG-Model' ) {
        # Summarize the beta coefficients from the posterior
        ms <- posterior::summarise_draws(obj$beta) |>
          subset(select = c('variable', 'mean'))
        # Reorder
        ms <- ms[order(ms$mean,decreasing = TRUE),] # Sort

        message(paste0(
          'Trained ',self$name,' (',self$show(),')',
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
      } else if( self$get_name() == 'GLMNET-Model') {
        # Summarise coefficients within 1 standard deviation
        ms <- tidy_glmnet_summary(obj)

        message(paste0(
          'Trained ',self$name,' (',self$show(),')',
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

      } else if( self$get_name() == 'GLM-Model' ) {
        # Summarise coefficients within 1 standard deviation
        ms <- tidy_glm_summary(obj)

        message(paste0(
          'Trained ',self$name,' (',self$show(),')',
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

      } else if( self$get_name() == 'SCAMPR-Model' ) {
        # Summarise coefficients within 1 standard deviation
        ms <- obj$fixed.effects |>
          as.data.frame() |>
          tibble::rownames_to_column(var = "variable")
        # Remove intercept
        int <- grep("Intercept",ms$variable,ignore.case = TRUE)
        if(length(int)>0) ms <- ms[-int,]

        # Rename the estimate and std.error column
        ms <- ms |> dplyr::rename(mean = "Estimate", se = "Std. Error")

        message(paste0(
          'Trained ',self$name,' (',self$show(),')',
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

    #' @description
    #' Show the name of the Model.
    #' @return A [`character`] of the run name.
    show = function() {
      assertthat::assert_that(is.list(self$model))
      return( self$model$runname )
    },

    #' @description
    #' Plots the prediction if found.
    #' @param what [`character`] with the specific layer to be plotted.
    #' @return A graphical representation of the prediction
    plot = function(what = 'mean'){
      if( length( self$fits ) != 0 && !is.null( self$fits$prediction ) ){
        pred <- self$get_data('prediction')
        assertthat::assert_that(is.Raster(pred))
        # Check if median is requested but not present, change to q50
        if(what == "median" && !(what %in% names(pred))) { what <- "q50" }

        # Match argument
        what <- match.arg(what, names(pred), several.ok = FALSE)
        assertthat::assert_that( what %in% names(pred),
                                 msg = paste0('Prediction type not found. Available: ', paste0(names(pred),collapse = '|')))
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

    #' @description
    #' Plots the thresholded prediction if found.
    #' @param what [`character`] or [`numeric`] for the layer to be plotted.
    #' @return A graphical representation of the thresholded prediction if found.
    plot_threshold = function(what = 1){
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
        invisible(self)
      }
    },

    #' @description
    #' Show model run time if settings exist
    #' @return A [`numeric`] estimate of the duration it took to fit the models.
    show_duration = function(){
      if(!is.Waiver(self$settings)) self$settings$duration()
    },

    #' @description
    #' Get effects or importance tables from model
    #' @param obj A [`character`] of which object to return.
    #' @return A [`data.frame`] summarizing the model, usually its coefficient.
    summary = function(obj = 'fit_best'){
      assertthat::assert_that(
        is.character(obj)
      )
      # Distinguishing between model types
      if(self$get_name() ==  'GDB-Model'){
        clean_mboost_summary( self$get_data(obj) )
      } else if( self$get_name() == 'INLA-Model' || self$get_name() == 'INLABRU-Model'){
        tidy_inla_summary(self$get_data(obj))
      } else if( self$get_name() == 'BART-Model'){
        # Number of times each variable is used by a tree split
        # Tends to become less informative with higher numbers of splits
        varimp.bart(self$get_data(obj)) |> tibble::remove_rownames()
      } else if( self$get_name() == 'STAN-Model'){
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
      } else if( self$get_name() == 'BREG-Model'){
        posterior::summarise_draws(self$get_data(obj)$beta)
      } else if( self$get_name() == "XGBOOST-Model"){
        xgboost::xgb.importance(model = self$get_data(obj))
      } else if( self$get_name() == 'GLMNET-Model'){
        tidy_glmnet_summary(self$get_data(obj))
      } else if( self$get_name() == 'GLM-Model'){
        tidy_glm_summary(self$get_data(obj))
      } else if( self$get_name() == 'SCAMPR-Model'){
        # Summarise coefficients within 1 standard deviation
        ms <- self$get_data(obj)$fixed.effects |>
          as.data.frame() |>
          tibble::rownames_to_column(var = "variable")
        # Remove intercept
        int <- grep("Intercept",ms$variable,ignore.case = TRUE)
        if(length(int)>0) ms <- ms[-int,]
        # Rename the estimate and std.error column
        ms |> dplyr::rename(mean = "Estimate", se = "Std. Error")
      }
    },

    #' @description
    #' Generic plotting function for effect plots
    #' @param x A [`character`] for the object in question.
    #' @param what A [`character`] for the type of coefficients.
    #' @param ... Any other options.
    #' @return A graphical representation of the coefficents.
    effects = function(x = 'fit_best', what = 'fixed', ...){
      assertthat::assert_that(is.character(what))
      # Get model
      obj <- self$get_data(x)
      if( self$get_name() == 'GDB-Model'){
        # How many effects
        n <- length( stats::coef( obj ))
        # Use the base plotting
        par.ori <- graphics::par(no.readonly = TRUE)
        graphics::par(mfrow = c(ceiling(n/3),3))
        mboost:::plot.mboost(x = obj, type = 'b',cex.axis=1.5, cex.lab=1.5)
        graphics::par(par.ori)#dev.off()
      } else if( self$get_name() == 'INLA-Model') {
        plot_inla_marginals(obj, what = what)
      } else if( self$get_name() == 'GLMNET-Model') {
        if(what == "fixed"){
          ms <- tidy_glm_summary(obj)
          graphics::dotchart(ms$mean,
                             labels = ms$variable,
                             frame.plot = FALSE,
                             color = "grey20")
        } else{ plot(obj) }
      } else if( self$get_name() == 'GLM-Model') {
        if(what == "fixed"){
          glmnet:::plot.glmnet(obj$glmnet.fit, xvar = "lambda") # Deviance explained
        } else{ plot(obj) }
      } else if( self$get_name() == 'STAN-Model') {
        # Get true beta parameters
        ra <- grep("beta", names(obj),value = TRUE) # Get range
        rstan::stan_plot(obj, pars = ra)
      } else if( self$get_name() == 'INLABRU-Model') {
        # Use inlabru effect plot
        ggplot2::ggplot() +
          inlabru::gg(obj$summary.fixed, bar = TRUE)
      } else if( self$get_name() == 'BART-Model'){
        message('Calculating partial dependence plots')
        self$partial(obj, x.var = what, ...)
      } else if( self$get_name() == 'BREG-Model'){
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
      } else if( self$get_name() == "XGBOOST-Model"){
        # Check whether linear model was fitted, otherwise plot tree
        if( self$settings$get("only_linear") ){
          vi <- self$summary(x)
          xgboost::xgb.ggplot.importance(vi)
        } else {
          xgboost::xgb.plot.multi.trees(obj)
        }
      } else if( self$get_name() == "SCAMPR-Model"){
        dotchart(obj$fixed.effects[,1])
      } else {
        self$partial(obj, x.var = NULL)
      }
    },

    #' @description
    #' Get equation
    #' @return A [`formula`] of the inferred model.
    get_equation = function(){
      self$get_data("fit_best_equation")
    },

    #' @description
    #' Get specific fit from this Model
    #' @param x A [`character`] stating what should be returned.
    #' @return A [`SpatRaster`] object with the prediction.
    get_data = function(x = "prediction") {
      rr <- names(self$fits)
      if(!x %in% names(self$fits)){
        # Check if x is present in rr, if so print a message
        if(length(grep(x,rr))>0){
          if(getOption('ibis.setupmessages', default = TRUE)){
            myLog('[Estimation]','yellow','Output not found, but found: ', grep(x,rr,value = TRUE)[1])
          }
        }
        return(new_waiver())
      }
      return(self$fits[[x]])
    },

    #' @description
    #' Set new fit for this Model.
    #' @param x The name of the new fit.
    #' @param value The [`SpatRaster`] layer to be inserted.
    #' @return This object.
    set_data = function(x, value) {
      # Get projected value
      ff <- self$fits
      # Set the object
      ff[[x]] <- value
      self$fits <- ff
      return( self )
    },

    #' @description
    #' Get the threshold value if calculated
    #' @return A [`numeric`] threshold value.
    get_thresholdvalue = function(){
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

    #' @description
    #' Get threshold type and format if calculated.
    #' @return A vector with a [`character`] method and [`numeric`] threshold value.
    get_thresholdtype = function(){
      # Determines whether a threshold exists and plots it
      rl <- self$show_rasters()
      if(length(grep('threshold',rl))==0) return( new_waiver() )

      # Get the thresholded layer and return the respective attribute
      obj <- self$get_data( grep('threshold',rl,value = TRUE) )
      assertthat::assert_that(
        assertthat::has_attr(obj, "format"),
        assertthat::has_attr(obj, "method"))
      return(
        c("method" = attr(obj, "method"), "format" = attr(obj, "format"))
      )
    },

    #' @description
    #' List all rasters in object
    #' @return A [`vector`] with [`logical`] flags for the various objects.
    show_rasters = function(){
      rn <- names(self$fits)
      rn <- rn[ which( sapply(rn, function(x) is.Raster(self$get_data(x)) ) ) ]
      return(rn)
    },

    #' @description
    #' Get projection of the background.
    #' @return A geographic projection
    get_projection = function(){
      sf::st_crs(self$model$background)
    },

    #' @description
    #' Get the resolution of the projection
    #' @return [`numeric`] estimates of the distribution.
    get_resolution = function(){
      if(!is.Waiver(self$get_data())){
        terra::res( self$get_data() )
      } else {
        # Try to get it from the modelling object
        self$model$predictors_object$get_resolution()
      }
    },

    #' @description
    #' Remove calculated thresholds
    #' @return Invisible
    rm_threshold = function(){
      rl <- self$show_rasters()
      if(length(grep('threshold',rl))>0){
        for(val in grep('threshold',rl,value = TRUE)){
          self$fits[[val]] <- NULL
        }
      }
      invisible(self)
    },

    #' @description
    #' Calculate a suitability index for a given projection
    #' @details
    #' Methods can either be normalized by the minimum and maximum.
    #' Or the relative total using the sumof values.
    #' @param method The method used for normalization.
    #' @return Returns a [`SpatRaster`].
    calc_suitabilityindex = function(method = "normalize"){
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

    #' @description
    #' Get centroids of prediction layers
    #' @param patch A [`logical`] if centroid should be calculated weighted by values.
    #' @param layer [`character`] of the layer to use.
    #' @return Returns a [`sf`] object.
    get_centroid = function(patch = FALSE, layer = "mean"){
      assertthat::assert_that(
        is.logical(patch),
        is.character(layer),
        is.Raster(self$get_data()),
        msg = "This function only works if predictions have been calculated!"
      )

      # Check if threshold is present
      rl <- self$show_rasters()
      if(length( grep('threshold',rl,value = TRUE) )>0){
        # Threshold present
        obj <- self$get_data( grep('threshold',rl,value = TRUE) )
        assertthat::assert_that(length(grep(layer, names(obj),value = TRUE))>0)
        obj <- obj[[grep(layer, names(obj),value = TRUE)]]
        # Get format
        if(attr(obj, "format") == "binary"){
          # Calculate centroid
          cent <- raster_centroid(obj, patch)
        } else {
          # In this case we assume all values larger than 0 to be genuine patches
          # TODO: This could be implemented if there is need.
          obj[obj>0] <- 1
          cent <- raster_centroid(obj, patch)
        }
      } else {
        # Get non-thresholded layer
        obj <- self$get_data( grep('prediction',rl,value = TRUE) )
        assertthat::assert_that(layer %in% names(obj))
        obj <- obj[[layer]]
        # Calculate centroid, patch to FALSE as this is non-sensical here
        cent <- raster_centroid(obj, patch = FALSE)
      }
      return(cent)
    },

    #' @description
    #' Logical indication if the prediction was limited.
    #' @return A [`logical`] flag.
    has_limits = function(){
      # Check for settings
      settings <- self$settings
      if(!is.Waiver(settings)){
        return(
          settings$get('has_limits')
        )
      }
    },

    #' @description
    #' Has a offset been used?
    #' @return A [`logical`] flag.
    has_offset = function(){
      model <- self$model$offset
      if(!is.Waiver(model$offset)) return( TRUE )
      # Also check whether offset is somehow in the equation
      ind <- attr(stats::terms.formula(fit$get_equation()), "offset")
      if(!is.null(ind)) return( TRUE )
    },

    #' @description
    #' Convenience function to mask all input datasets.
    #' @param mask A \code{SpatRaster} or `sf` object.
    #' @param inverse A `logical` flag if the inverse should be masked instead.
    #' @param ... Any other parameters passed on to mask
    #' @return Invisible
    mask = function(mask, inverse = FALSE, ...){
      # Check whether prediction has been created
      prediction <- self$get_data()
      if(!is.Waiver(prediction)){
        # If mask is sf, rasterize
        if(inherits(mask, 'sf')){
          mask <- terra::rasterize(mask, prediction)
        }
        # Check that mask aligns
        if(!terra::compareGeom(prediction,mask)){
          mask <- terra::resample(mask, prediction, method = "near")
        }
        # Now mask and save
        prediction <- terra::mask(prediction, mask, inverse = inverse, ...)

        # Save data
        self$fits[["prediction"]] <- prediction

        # Do the same for any thresholds eventually found
        tr <- grep("threshold", self$show_rasters(), value = TRUE)
        if(length(tr)){
          m <- self$get_data(x = tr)
          m <- terra::mask(m, mask, inverse = inverse, ...)
          self$fits[[tr]] <- m
        }
        invisible(self)
      }
    },

    #' @description
    #' Save the prediction as output.
    #' @param fname An output filename as [`character`].
    #' @param type A format as [`character`]. Matched against a list of supported formats.
    #' @param dt The datatype used, such as float64
    #' @return Saved spatial prediction on drive.
    save = function(fname, type = 'gtif', dt = 'FLT4S'){
      assertthat::assert_that(
        is.character(fname),
        type %in% c('gtif','gtiff','tif','nc','ncdf'),
        'fits' %in% names(self),
        dt %in% c('LOG1S','INT1S','INT1U','INT2S','INT2U','INT4S','INT4U','FLT4S','FLT8S')
      )
      type <- tolower(type)
      if(type %in% c("gtif", "gtiff", "tif")){
        fname <- paste0(tools::file_path_sans_ext(fname), ".tif")
      }

      # Get raster file in fitted object
      cl <- sapply(self$fits, class)
      if(length( grep('SpatRaster', cl,ignore.case = T) )==0){
        # Security check in case for some reason there are no predictions
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Output]','red','No predictions found?')
        return(NULL)
      }
      ras <- self$fits[grep('SpatRaster', cl,ignore.case = T)]
      assertthat::assert_that(length(ras)>0,
                              msg = "No prediction to save found.")

      # If is a list (multiple SpatRaster) -> Combine
      if(is.list(ras)) ras <- Reduce('c', ras)

      # Check that no-data value is not present in ras
      assertthat::assert_that(any(!terra::global(ras, "min", na.rm = TRUE)[,1] <= -9999),
                              msg = 'No data value -9999 is potentially in prediction!')

      if(file.exists(fname)) warning('Overwritting existing file...')
      if(type %in% c('gtif','gtiff','tif')){
        # Save as geotiff
        writeGeoTiff(ras, fname = fname, dt = dt)
      } else if(type %in% c('nc','ncdf')) {
        # Save as netcdf
        writeNetCDF(ras, fname = fname, varName = 'ibis.iSDM prediction',
                    varUnit = "Suitability",varLong = "Relative suitable habitat")
      }
      invisible(self)
    }

  ),

  # Any private entries
  private = list(
    finalize = function() {
    }
  ),
  # Don't lock objects so that engine-specific functions can be added
  lock_objects = FALSE,
  lock_class = FALSE
)

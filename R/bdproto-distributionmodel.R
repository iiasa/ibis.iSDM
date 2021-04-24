#' @include utils.R bdproto.R utils-inla.R
NULL

#' @export
if (!methods::isClass("DistributionModel")) methods::setOldClass("DistributionModel")
NULL

#' Prototype for the trained Model object
#'
#' All trained Models should inherit the options here
#'
#' @name DistributionModel-class
#'
#' @aliases DistributionModel
NULL

#' @export
DistributionModel <- bdproto(
  "DistributionModel",
  model = list(),
  fits = list(), # List of fits with data
  # Print message with summary of model
  print = function(self) {
    # TODO: Have a lot more information in here and to be prettified
    # FIXME: Have engine-specific code moved to engine
    if(inherits(self, 'INLA-Model')){
      if( length( self$fits ) != 0 ){
        # Get strongest effects
        ms <- subset(tidy_inla_summary(self$get_data('fit_best')),
                     select = c('variable', 'mean'))
        ms <- ms[order(ms$mean,decreasing = TRUE),] # Sort

        # TODO: See if this can be plotted and extracted and shown here
        # Provides the distance value (in the unit of the point coordinates) above
        # which spatial dependencies become negligible
        # model0.res<-inla.spde2.result(model0, 'spatial.field', spde, do.transf=TRUE)
        # model0.res$summary.log.range.nominal

        message(paste0(
          'Trained ',class(self)[1],' (',self$model$name,')',
          '\n  \033[2mStrongest summary effects:\033[22m',
          '\n     \033[34mPositive:\033[39m ', name_atomic(ms$variable[ms$mean>0]),
          '\n     \033[31mNegative:\033[39m ', name_atomic(ms$variable[ms$mean<0])
        ))
      }
    } else if( inherits(self, 'GDB-Model') ) {

        # Get Variable importance
        vi <- mboost::varimp(
          self$get_data('fit_best')
        )
        vi <- sort( vi[which(vi>0)],decreasing = TRUE )

        message(paste0(
          'Trained ',class(self)[1],' (',self$model$name,')',
          '\n  \033[2mStrongest effects:\033[22m',
          '\n     ', name_atomic(names(vi))
        ))
      } else {
      message(paste0(
        'Trained distribution model (',self$model$name,')',
        '\n     No fitted model found!'
      ))
    }
  },
  # Show the name of the Model
  show = function(self) {
    self$model$name
  },
  # Plot the prediction
  plot = function(self, what = 'mean',...){
    if( length( self$fits ) != 0 && !is.null( self$fits$prediction ) ){
      pred <- self$get_data('prediction')
      assertthat::assert_that(
        inherits(pred,'Raster'), what %in% names(pred)
      )
      # Custom colour column
      cols <- colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'))(100)

      raster::plot(pred[[what]],
           main = paste0(what, ' prediction (',what,')'),
           box = FALSE,
           axes = TRUE,
           colNA = NA, col = cols,...)

    } else {
      message(
        paste0('No model predictions found.')
      )
    }
  },
  # Get effect tables from model
  summary = function(self, x = 'fit_best'){
    # Distinguishing between model types
    if(inherits(self, 'GDB-Model')){
      mboost:::summary.mboost(mod1$get_data(x))
    } else if(inherits(self, 'INLA-Model')){
      tidy_inla_summary(self$get_data(x))
    }
  },
  # Generic plotting function for partial effects
  effects = function(self, x = 'fit_best', what = 'fixed'){
    if(inherits(self, 'GDB-Model')){
      # How many effects
      n <- length( coef( self$get_data(x) ))
      # Use the base plotting
      par.ori <- par(no.readonly = TRUE)
      par(mfrow = c(ceiling(n/3),3))

      mboost:::plot.mboost(x = self$get_data(x),
                           type = 'b',cex.axis=1.5, cex.lab=1.5)

      par(par.ori)#dev.off()
    } else if(inherits(self, 'INLA-Model')) {
      plot_inla_marginals(self$get_data(x),what = 'fixed')
    }
  },
  # Get specific fit from this Model
  get_data = function(self, x) {
    if (!x %in% names(self$fits))
      return(new_waiver())
    return(self$fits[[x]])
  },
  # Set fit for this Model
  set_data = function(self, x, value) {
    self$fits[[x]] <- value
    invisible()
  }
)

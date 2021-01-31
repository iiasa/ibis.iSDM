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
  name = character(0),
  id = new_waiver(),
  fits = list(), # List of fits with data
  # Print message with summary of model
  print = function(self) {
    # TODO: Have a lot more information in here and to be prettified
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
        'Trained distribution model (',self$name,')',
        '\n  \033[2mStrongest summary effects:\033[22m',
        '\n     \033[34mPositive:\033[39m ', name_atomic(ms$variable[ms$mean>0]),
        '\n     \033[31mNegative:\033[39m ', name_atomic(ms$variable[ms$mean<0])
      ))

    } else {
      message(paste0(
        'Trained distribution model (',self$name,')',
        '\n     No fitted model found!'
      ))
    }
  },
  # Show the name of the Model
  show = function(self) {
    self$name
  },
  # Plot the prediction
  plot = function(self, what = 'mean',...){
    if( length( self$fits ) != 0 ){
      pred <- self$get_data('prediction')
      assertthat::assert_that(
        inherits(pred,'Raster'), what %in% names(pred)
      )
      # Magma column
      cols <- c("#0D0887FF","#47039FFF","#7301A8FF","#9C179EFF","#BD3786FF","#D8576BFF",
                "#ED7953FF","#FA9E3BFF","#FDC926FF","#F0F921FF")
      raster::plot(pred[[what]],
           main = paste0(what, ' prediction (',self$name,')'),
           box = FALSE,
           axes = TRUE,
           colNA = NA, col = cols,...)

    } else {
      message(
        paste0('No trained models found.')
      )
    }
  },
  # Plot SPDE if existing
  plot_spde = function(self,...){
    if( length( self$fits$fit_best$size.spde2.blc ) == 1)
      {
      # Get spatial projections from model
      # FIXME: Potentially make the plotting of this more flexible
      gproj <- INLA::inla.mesh.projector(self$get_data('mesh'),  dims = c(300, 300))
      g.mean <- INLA::inla.mesh.project(gproj,
                                        self$get_data('fit_pred')$summary.random$spatial.field$mean)
      g.sd <- INLA::inla.mesh.project(gproj, self$get_data('fit_pred')$summary.random$spatial.field$sd)

      # Plot
      par(mfrow=c(1,2))
      raster::image(g.mean,main = 'mean spatial effect')
      raster::image(g.sd,main = 'sd spatial effect')


      } else {
        message('No spatial covariance in model specified.')
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

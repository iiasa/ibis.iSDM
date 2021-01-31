#' Function for creating a joint fitted and prediction stack
#' @param stk_resp A stack object
#' @param cov The covariate data stack
#' @param mesh The background projection mesh
#' @param type Name to use
#' @param spde An spde field if specified
#' @noRd

inla_make_prediction_stack <- function(stk_resp, cov, mesh, type, spde = NULL){
  # Security checks
  assertthat::assert_that(
    inherits(stk_resp, 'inla.data.stack'),
    inherits(mesh,'inla.mesh'),
    is.character(type),
    is.data.frame(cov),
    is.null(spde)  || inherits(spde,'inla.spde')
  )
  # need to create a new A projection matrix for the prediction data
  mat_pred <- INLA::inla.spde.make.A(mesh,
                             loc = as.matrix(cov[,1:2]) )

  # Single response
  ll_pred <- list()
  ll_pred[[ stk_resp$data$names[[1]] ]] <- NA # Set to NA to predict for fitted area
  # Effects
  ll_effects <- list()
  # Note, order adding this is important apparently...
  ll_effects[['intercept']] <- list(intercept = rep(1,mesh$n) ) # FIXME: Potential source for bug. Think name of intersects need to differ if multiple models specified
  ll_effects[['predictors']] <- cov
  if(!is.null(spde)) ll_effects[['spatial.field']] <- list(Bnodes = 1:spde$n.spde)
  # Define A
  if(!is.null(spde)) {
      A = list(mat_pred, 1, mat_pred)
  } else {
      A = list(mat_pred, 1)
  }

  # Create stack depending on the number of variables in response
  if( length(stk_resp$data$names[[1]]) > 1) {
    stop('TBD')
    ys <- cbind(rep(NA, nrow(pred.grid)), rep(NA, nrow(pred.grid)))
    stack.pred.response <- inla.stack(data=list(y=ys),
                                      effects = list(list(data.frame(interceptA=rep(1,np))), env = pred.grid$cov, list(uns_field=1:spde$n.spde)),
                                      A = A,
                                      tag = paste0('pred_', type))
  } else if("Ntrials" %in% stk_resp$data$names) {
    stop('TBD')
    stack.pred.response <- inla.stack(data=list(y=NA, Ntrials = rep(1,np)),
                                      effects = list(list(data.frame(interceptA=rep(1,np))), env = pred.grid$cov, list(Bnodes=1:spde$n.spde)),
                                      A = A,
                                      tag = paste0('pred_', type) )
  } else {
    stk_pred <-
      INLA::inla.stack(
        data =  ll_pred,              # Response
        A = A,                        # Predictor projection matrix
        effects = ll_effects,         # Effects matrix
        tag = paste0('pred_',type) # New description tag
      )
  }

  # TODO: Insert some security checks that stk_resp and stK_pred are equal in variable names
  # Final security checks
  # assertthat::assert_that(
  #   names(stk_resp$effects$names) %in% names(stk_pred$effects$names),
  #   msg = 'Response stak and prediction stack mismatch.'
  # )
  return(stk_pred)
}

#' Tidy up summary information from a INLA model
#'
#' @param m A trained INLA model object
#' @param what A [`character`]
#' @param ... Other options to based on
#' @noRd

#TODO: Lot more to add here, including options on what to extract
tidy_inla_summary <- function(m, what = 'fixed',...){
  assertthat::assert_that(
    inherits(m,'inla'),
    is.character(what),length(what)==1,
    what %in% c('fixed','fitted','random','spde2')
  )

  w1 <- grep('summary',names(m),value = TRUE) # Grep summary objects
  w2 <- grep(what, w1,value = TRUE) # Grep specific summary
  assertthat::assert_that(length(w2)==1)

  # Format the output
  m[[w2]] %>%
    tibble::rownames_to_column('variable') %>%
    tibble::as_tibble()
}

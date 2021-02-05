#' @include utils.R waiver.R bdproto.R
NULL

#' @export
if (!methods::isClass("PredictorDataset")) methods::setOldClass("PredictorDataset")
NULL

#' PredictorDataset prototype description
#'
#' @name PredictorDataset-class
#' @aliases PredictorDataset
NULL

#' @export
PredictorDataset <- bdproto(
  "PredictorDataset",
  id           = character(0),
  data         = new_waiver(),
  # Printing function
  # FIXME: Prettify. Print summary values for each layer
  print = function(self){
    message(paste0(self$name(),':',
                   '\n Name(s):  ',name_atomic(self$get_names(), "predictors")
    ))
  },
  # Return name
  name = function(self){
    'Predictor dataset'
  },
  # Get Id
  id = function(self){
    self$id
  },
  # Get names of data
  get_names = function(self){
    names(self$get_data())
  },
  # Get data
  get_data = function(self, df = FALSE, na.rm = TRUE, ...){
    if(df) {
      as.data.frame(self$data, xy = TRUE,na.rm = na.rm, ...)
    } else self$data
  },
  # Add a new Predictor dataset to this collection
  set_data = function(self, x, value){
    assertthat::assert_that(assertthat::is.string(x),
                            inherits(value, "Raster"))
    self$data <- addLayer(self$get_data(),value)
    invisible()
  },
  # Remove a specific Predictor by name
  rm_data = function(self, x) {
    assertthat::assert_that(is.vector(x) || is.character(x),
                            all(x %in% names(self$get_data()))
                            )
    # Match indices
    ind <- match(x, self$get_names())
    # Overwrite predictor dataset
    self$data <- raster::dropLayer(self$get_data(), ind)
    invisible()
  },
  # Print input messages
  show = function(self) {
    self$print()
  },
  # Collect info statistics with optional decimals
  get_summary = function(self, digits = 2) {
    # Maybe make a little bit prettier
    round(
      summary( self$get_data()), digits = digits
         )
  },
  # Number of Predictors in object
  length = function(self) {
    if(inherits(self$get_data(),'Raster'))
      nlayers(self$get_data())
    else
      ncol(self$get_data)
  }
)

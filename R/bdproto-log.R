#' @include utils.R bdproto.R
NULL

#' @export
if (!methods::isClass("Log")) methods::setOldClass("Log")
NULL

#' Log prototype.
#'
#' @description
#' Basic proto object for Log, any Log inherit from here
#'
#' @family bdproto
#' @keywords bdproto
#'
#' @name Log-class
NULL

#' @rdname Log-class
#' @export
Log <- bdproto(
  "Log",
  filename = character(0),
  output = new_waiver(),
  # Print message with filename
  print = function(self) {
    message( paste0('Log: ', paste0('\033[45m', basename(self$filename), '\033[49m' ) ) )
  },
  # Opens the connection to the output filename
  open = function(self){
    assertthat::assert_that(
      assertthat::has_extension(self$filename,'txt')
    )
    try({output = file( self$filename, open = 'wt', method = 'default') },silent = TRUE)

    if(class(output)[1]=='try-error') { stop('Log file could not be opened!')}

    sink(output, append = TRUE, type = "output") # Writing console output to log file
    sink(output, append = TRUE, type = "message") # Writing messages to the log file

    self$output <- output
    invisible()
  },
  # Closes the connection to the output file
  close = function(self) {
    assertthat::assert_that(!is.Waiver(self$output))
    # Close connection to log file
    closeAllConnections()
    # Save waiver again
    self$output <- new_waiver()
    invisible()
  },
  # Get output filename
  get_filename = function(self){
    return( paste0( basename(self$filename)  ) )
  },
  # Set output filename
  set_filename = function(self, value){
    assertthat::assert_that(is.character(value),
                            assertthat::has_extension(value,'txt')
                            )
    self$filename <- value
    invisible()
  },
  # Delete log file
  delete = function(self){
    file.remove(self$filename)
    invisible()
  },
  # Open with system viewer
  open_system = function(self){
    assertthat::assert_that(file.exists(self$filename))
    shell.exec(self$get_filename())
    invisible()
  }
)

if (!methods::isClass("Log")) methods::setOldClass("Log")

#' Log prototype.
#'
#' @description
#' Basic [`R6`] object for Log, any Log inherit from here
#'
#' @keywords classes
#'
#' @name Log-class
NULL

#' @rdname Log-class
#' @export
Log <- R6::R6Class(
  "Log",
  public = list(
    #' @field filename A [`character`] of where the log is to be stored.
    #' @field output The log content.
    filename = character(0),
    output = new_waiver(),

    #' @description
    #' Initializes the object and specifies some default parameters.
    #' @param filename A [`character`] of where the log is to be stored.
    #' @param output The log content.
    #' @return NULL
    initialize = function(filename, output){
      assertthat::assert_that(
        is.character(filename)
      )
      self$filename <- filename
      self$output <- output
    },

    #' @description
    #' Print message with filename
    #' @return A message on screen
    print = function() {
      message( paste0('Log: ', paste0('\033[45m', basename(self$filename), '\033[49m' ) ) )
    },

    #' @description
    #' Opens the connection to the output filename.
    #' @param type A [`character`] vector of the output types.
    #' @return Invisible TRUE
    open = function(type = c("output", "message")){
      assertthat::assert_that(
        assertthat::has_extension(self$filename,'txt'),
        length(type)>0,
        all( sapply(type, is.character) )
      )
      try({output = base::file( self$filename, open = 'wt', method = 'default') },
          silent = TRUE)

      if(class(output)[1]=='try-error') { stop('Log file could not be opened!')}

      # For each type append
      # Writing console or message output to log file
      for(tt in type){
        base::sink(output, append = TRUE, type = tt)
      }

      self$output <- output
      invisible()
    },

    #' @description
    #' Closes the connection to the output file
    #' @return Invisible TRUE
    close = function() {
      assertthat::assert_that(!is.Waiver(self$output))
      # Close connection to log file
      base::closeAllConnections()
      # Save waiver again
      self$output <- new_waiver()
      invisible()
    },

    #' @description
    #' Get output filename
    #' @return A [`character`] with the filename
    get_filename = function(){
      return( paste0( basename(self$filename)  ) )
    },

    #' @description
    #' Set a new output filename
    #' @param value A [`character`] with the new filename.
    #' @return Invisible TRUE
    set_filename = function(value){
      assertthat::assert_that(is.character(value),
                              assertthat::has_extension(value,'txt')
      )
      self$filename <- value
      invisible()
    },

    #' @description
    #' Delete log file
    #' @return Invisible TRUE
    delete = function(){
      if(file.exists(self$filename)) file.remove(self$filename)
      invisible()
    },

    #' @description
    #' Open log with system viewer
    #' @return Invisible TRUE
    open_system = function(){
      assertthat::assert_that(file.exists(self$filename))
      shell.exec(self$get_filename())
      invisible()
    }
  ),

  # Any private entries
  private = list(
    finalize = function() {
    }
  )
)

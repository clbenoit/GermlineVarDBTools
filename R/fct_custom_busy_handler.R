#' buildDB_seqone
#' @description
#' A short description...
#' @param name description
#' @return The return value, if any, from executing the function.
#'
#' @noRd
#' @export
custom_busy_handler <- function(retries = 60) {
  if (retries < 100) {  # Retry up to 100 times
    Sys.sleep(0.1)  # Sleep for 100 milliseconds between retries
    return(100)  # Indicate 100 milliseconds delay before retrying
  }
  return(0)  # Stop retrying
}
  
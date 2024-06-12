#' buildDB_seqone
#' @description
#' A short description...
#' @param name description
#' @importFrom dplyr %>% select filter mutate inner_join group_by
#' @importFrom RSQLite SQLite
#' @importFrom DBI dbGetQuery dbExistsTable dbWriteTable dbDisconnect dbConnect
#' @return The return value, if any, from executing the function.
#'
#' @noRd
#' @export
initcomment <- function(prefix = NULL, db_path = NULL) {

  db_name <- file.path(db_path,paste0(prefix, ".db"))
  ## add samples ##
  con <- dbConnect(SQLite(), db_name)
  if(!(DBI::dbExistsTable(conn = con, name = "variant_comments"))){
    new_comment <- data.frame(com_id = 1, variant_id = 'NULL', comment = 'NULL' , user = 'NULL',date = 'NULL')
    DBI::dbWriteTable(conn = con, name="variant_comments", new_comment, overwrite = TRUE)
  }
  
 if(!(DBI::dbExistsTable(con,"sample_comments"))){
   new_comment <- data.frame(com_id = 1, sample = 'NULL', comment = 'NULL', user = 'NULL',date = 'NULL')
   DBI::dbWriteTable(conn = con, name="sample_comments",value = new_comment, overwrite = TRUE)
  }
  
}

#' QCtodb_seqone
#' @description
#' A short description...
#' @param name description
#' @importFrom dplyr %>% select filter mutate inner_join group_by rename
#' @importFrom RSQLite SQLite
#' @importFrom DBI dbGetQuery dbExistsTable dbWriteTable dbDisconnect dbConnect
#' @return The return value, if any, from executing the function.
#'
#' @noRd
#' @export
QCtodb_seqone <- function(prefix = NULL,
                          qcfile_path = NULL,
                          db_path = NULL) {
  
  ## add samples ##
  db_name <- file.path(db_path,paste0(prefix, ".db"))
  con <- dbConnect(SQLite(), db_name)
  
  qc_new <- read.table(qcfile_path,header = TRUE, sep = "\t", check.names = FALSE)  %>% rename(sample = Sample)# %>% select(c())
  
  if(DBI::dbExistsTable(conn = con, name="QC")){
    qc <- DBI::dbReadTable(conn = con, name="QC", check.names = FALSE)
    qc_new <- qc_new[!(qc_new$sample %in% qc$sample),] 
  } else { qc <- qc_new }

  DBI::dbWriteTable(conn = con,
                   name = "QC", 
                   value = qc_new,
                   append = TRUE)
  
  DBI::dbDisconnect(con)
  
}

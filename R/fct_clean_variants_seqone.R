#' compute_frequencies
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#' @importFrom dplyr %>% select filter mutate inner_join group_by
#' @importFrom RSQLite SQLite
#' @importFrom DBI dbGetQuery dbExistsTable dbWriteTable dbDisconnect dbConnect
#' @noRd
#' @export
clean_variants_seqone <- function(db_path = NULL,
                    attribute = NULL, prefix = NULL){

 
  db_name <- file.path(db_path,paste0(prefix, ".db"))
  con <- DBI::dbConnect(RSQLite::SQLite(), db_name)

  dbSendQuery(con,"DELETE FROM variant_geno WHERE qa/dp <15;")
  try({dbSendQuery(con,"DROP INDEX idx_variant_geno;")})
  dbSendQuery(con,"CREATE INDEX idx_variant_geno ON variant_geno (variant_id);")
  
  dbSendQuery(con,"DELETE FROM variant_impact WHERE variant_id NOT IN (SELECT variant_id FROM variant_geno);")
  try({dbSendQuery(con,"DROP INDEX idx_variant_impact;")})
  dbSendQuery(con,"CREATE INDEX idx_variant_impact ON variant_impact (variant_id);")
  
  dbSendQuery(con,"DELETE FROM variant_info WHERE variant_id NOT IN (SELECT variant_id FROM variant_geno);")
  try({dbSendQuery(con,"DROP INDEX idx_variant_info;")})
  dbSendQuery(con,"CREATE INDEX idx_variant_info ON variant_info (variant_id);")
  
  dbDisconnect(con)
  
  

}



#' compute_frequencies
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#' @importFrom RSQLite SQLite
#' @importFrom DBI dbSendQuery dbDisconnect dbConnect dbDisconnect
#' @noRd
#' @export
rebuild_indexes <- function(db_path = NULL,
                    attribute = NULL, prefix = NULL){
  
  db_name <- file.path(db_path,paste0(prefix, ".db"))
  con <- DBI::dbConnect(RSQLite::SQLite(), db_name)

  try({dbSendQuery(con,"DROP INDEX idx_variant_geno;")})
  dbSendQuery(con,"CREATE INDEX idx_variant_geno ON variant_geno (variant_id);")
  
  try({dbSendQuery(con,"DROP INDEX idx_variant_impact;")})
  dbSendQuery(con,"CREATE INDEX idx_variant_impact ON variant_impact (variant_id);")
  
  try({dbSendQuery(con,"DROP INDEX idx_variant_info;")})
  dbSendQuery(con,"CREATE INDEX idx_variant_info ON variant_info (variant_id);")
  
  dbDisconnect(con)

}

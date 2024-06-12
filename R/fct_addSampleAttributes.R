#' vcftodb
#' @description A fct function
#' @importFrom RSQLite SQLite
#' @importFrom DBI dbSendQuery dbDisconnect dbConnect dbDisconnect dbReadTable
#' @importFrom dplyr left_join
#' @noRd
#' @export
addSampleAttributes <- function(prefix = NULL, db_path = NULL, matches_file = NULL) {
  
  db_name <- file.path(db_path,paste0(prefix, ".db"))
  con <- dbConnect(SQLite(), db_name)
  matchs <- read.table(file  = matches_file, header = TRUE, sep = "\t")
  samples <- DBI::dbReadTable(con,"samples") %>% dplyr::select(-one_of(colnames(matchs)[colnames(matchs) != "sample"]))
  samples <- left_join(samples,matchs, by = "sample")
  
  DBI::dbWriteTable(conn = con,
                    name="samples",
                    value= samples,
                    overwrite = TRUE)
  
  new_attributes <- colnames(matchs)[!(colnames(matchs) %in% colnames(samples))]
  for (attribute in new_attributes){ # compute frequencies on DB for all new added attributes
    GermlineVarDB::compute_frequency(db_path = db_path, prefix = prefix ,attribute = attribute)
  }
}

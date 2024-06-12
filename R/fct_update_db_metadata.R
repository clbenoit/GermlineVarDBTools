#' compute_frequencies
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#' @importFrom digest digest
#' @noRd
#' @export
update_db_metadata <- function(db_path = NULL,
                    prefix = NULL){

  db_name <- file.path(db_path,paste0(prefix, ".db"))

  con <- dbConnect(SQLite(), db_name)
  samples <- DBI::dbReadTable(con,"samples")
  variant_geno <- DBI::dbReadTable(con,"variant_geno") #%>% # Besoin de tout load ici ? C'est juste pour les value de box et input
  variant_info <- DBI::dbReadTable(con,"variant_info")
  variants_list = unique(variant_info$variant_id)

  metadata = data.frame(dp_max = max(variant_geno$dp),
                        dp_min = min(variant_geno$dp),
                        af_max = max(as.numeric(variant_geno$af),na.rm = TRUE),
                        af_min = min(variant_geno$af,na.rm = TRUE),
                        qual_max = max(variant_geno$qa,na.rm = TRUE),
                        qual_min = min(variant_geno$qa,na.rm = TRUE),
                        nb_variants = as.character(length(variants_list)),
                        nb_samples = as.character(length(unique(samples$sample)))
                        )
  
  metadata$hash <- digest(metadata, algo = "md5")

  DBI::dbWriteTable(conn = con,
                    name="db_metadata",
                    value= metadata,
                    overwrite = TRUE)

}

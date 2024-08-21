#' compute_frequencies
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#' @importFrom dplyr full_join
#'
#' @noRd
#' @export
compute_sample_frequency <- function(db_path = NULL,
                    attribute = NULL, prefix = NULL){

  db_name <- file.path(db_path,paste0(prefix, ".db"))
  con <- dbConnect(SQLite(), db_name)
  
  variant_geno <- DBI::dbReadTable(con,"variant_geno")
  if(is.null(attribute)){
    frequencies <-aggregate(variant_geno$gt, by=list(Category=variant_geno$variant_id), FUN=sum)
    frequencies$x <- (frequencies$x/length(unique(variant_geno$sample)))
  } else {
    samples <- DBI::dbReadTable(con,"samples")
    samples <- samples %>% filter(sample %in% unique(variant_geno$sample))
    frequencies <- data.frame("variant_id" = 'a')
    columns <- DBI::dbGetQuery(conn = con, "PRAGMA table_info('frequencies');")$name    
    for (value in unique(samples[,attribute])){
      col <- columns[grepl(paste0(value,"_freq_total"),columns)]
      # print(paste0("ALTER TABLE frequencies DROP COLUMN ",col,";"))
      # tryCatch({
      #   DBI::dbSendQuery(conn = con, paste0("ALTER TABLE frequencies DROP COLUMN '",col,"';"))
      # }, error = function(e){ print(e) })
      
      samples_fil <- samples %>% filter(get(attribute) == value)
      rm(samples)
      variant_geno_filtered <- variant_geno %>% 
        select(c("variant_id","sample","gt")) %>%
        inner_join(samples_fil,
                   by = "sample", keep = FALSE)
      rm(variant_geno)
      variant_geno_filtered$gt <- 1
      n_samples <- length(unique(samples_fil$sample))
      frequencies_col <-aggregate(variant_geno_filtered$gt, by=list(variant_id=variant_geno_filtered$variant_id), FUN=sum)
      frequencies_col[,paste0(value,"_freq","_total_",n_samples)] <- (frequencies_col$x/n_samples) 
      
      frequencies_col <- frequencies_col %>% select(-x) 
      frequencies <- frequencies_col %>% full_join(frequencies, by = c("variant_id"), keep = FALSE)
    }
    frequencies[is.na(frequencies)] = 0
    if(!DBI::dbExistsTable(conn = con,name="frequencies")){
    DBI::dbWriteTable(conn = con, name="frequencies", value = frequencies, overwrite = TRUE) 
    } else {
    table <- DBI::dbReadTable(conn = con,  name="frequencies")
    table <- subset(table, select = which(!(names(table) %in% col)))
    frequencies <- left_join(frequencies,table, by = "variant_id", keep = FALSE)
    DBI::dbWriteTable(conn = con, name="frequencies", value = frequencies, overwrite = TRUE) 
    }
  }
  dbDisconnect(con)
}

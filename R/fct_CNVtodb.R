#' CNVtodbfinal_seqone
#' @description
#' A short description...
#' @param name description
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom RSQLite SQLite
#' @importFrom GenomicFeatures exonicParts exons
#' @importFrom DBI dbConnect dbExistsTable dbReadTable
#' @importFrom AnnotationDbi loadDb
#' @importFrom dplyr %>% select filter mutate inner_join group_by bind_rows summarise
#' @return The return value, if any, from executing the function.
#' @noRd
#' @export
CNVtodb_seqone <- function(prefix = NULL,
                          cnvfile_path = NULL,
                          db_path = NULL) {
  
  db_name <- file.path(db_path,paste0(prefix, ".db"))
  print(db_name)
  ## add samples ##
  con <- dbConnect(SQLite(), db_name)
  gtf <- loadDb(system.file("extdata","annotations/gencode.v43lift37.annotation_canonical", package = "GermlineVarDB"))
  exons <- GenomicFeatures::exonicParts(gtf)

  cnv_new_geno_raw <- read.table(cnvfile_path,header = TRUE, sep = "\t") %>% dplyr::select(c("sample_name","interval","padj",
                                                                                  "chr", "ratio", "gene", "refseq_accession_number",
                                                                                  "start","end", #"strand",
                                                                                  )) %>% mutate(chr = paste0("chr",chr))
  
  cnv_new_geno_granges <- GenomicRanges::makeGRangesFromDataFrame(cnv_new_geno_raw, ignore.strand = TRUE,keep.extra.columns = TRUE,
                                                                  start.field= "start", end.field = "end", seqnames.field = c("chr"))
  
  #### % Exon covered by of seqOneresult  ####
  hits <- findOverlaps(query = exons, subject = cnv_new_geno_granges, ignore.strand=TRUE, minoverlap = 2, type = "any")
  overlaps <- pintersect(exons[queryHits(hits)], cnv_new_geno_granges[subjectHits(hits)],ignore.strand=TRUE)
  percentOverlap <- width(overlaps) / width(cnv_new_geno_granges[subjectHits(hits)])
  hits <- hits[percentOverlap > 0.01]
  probes_covered_by_exon <- cnv_new_geno_granges[subjectHits(hits)]
  
  nrow(mcols(cnv_new_geno_granges))
  nrow(mcols(exons))
  nrow(mcols(probes_covered_by_exon))
  (nrow(mcols(probes_covered_by_exon))/nrow(cnv_new_geno_raw))*100
  
  ########### CEST DANS LETAPE SUIVANTE QUE LES LISTES SONT CREES, EVITABLE ? #############
  mcols(probes_covered_by_exon)$exon_rank <- mcols(exons)$exon_rank[queryHits(hits)]
  mcols(probes_covered_by_exon)$gene_id <- mcols(exons)$gene_id[queryHits(hits)]
  
  cnv_new_geno <- as.data.frame(probes_covered_by_exon) %>% dplyr::select(c("start","end","ratio","sample_name","gene","exon_rank","interval")) %>%
                     unnest(cols = c(exon_rank),.drop = FALSE) %>%
                     mutate(concat = paste(sample_name,gene,exon_rank, sep = "_"))  %>% group_by(concat) %>% 
    summarise(mean_ratio = as.numeric(mean(ratio)))  %>%
     mutate(CNV_status = case_when(mean_ratio <= 0.7   ~ "DEL",
                      mean_ratio >= 1.3   ~ "DUP",
                      mean_ratio >= 2   ~ "AMP",
                      (0.8 <= mean_ratio & mean_ratio <= 1.2) ~ "Unknown"
                      )) %>% 
    mutate(sample_name = stringr::str_extract(concat,"[^_]+")) %>%
    mutate(exon_rank = gsub("_","",stringr::str_extract(concat,"_[0-9]+"))) %>%
    mutate(gene = gsub("_","",stringr::str_extract(concat,"_.*_"))) %>%
    mutate(variant_id = paste(gene,exon_rank,CNV_status,sep = "_"))

  print("saving cnv genotypes")
  if(DBI::dbExistsTable(conn = con, name="cnv_geno")){
    cnv_geno <- as.data.frame(DBI::dbReadTable(conn = con,
                                 name="cnv_geno"))
    cnv_geno <- cnv_geno %>% filter(!(concat %in% cnv_new_geno$concat))
    cnv_new_geno <- as.data.frame(cnv_new_geno)
    cnv_new_geno <- bind_rows(cnv_geno, cnv_new_geno)
  }
  DBI::dbWriteTable(conn = con, value = cnv_new_geno, name="cnv_geno", append = TRUE)
  
  print("saving cnv raw genotypes")
  if(DBI::dbExistsTable(conn = con, name="cnv_geno_raw")){
    cnv_geno_raw <- DBI::dbReadTable(conn = con,
                                 name="cnv_geno_raw")
    cnv_geno_raw <- cnv_geno_raw %>% filter(!(sample_name %in% cnv_new_geno_raw$sample_name))
    cnv_new_geno_raw <- bind_rows(cnv_new_geno_raw,cnv_geno_raw)
  }
  DBI::dbWriteTable(conn = con, value = cnv_new_geno_raw, name="cnv_geno_raw", append = TRUE)
  
  DBI::dbDisconnect(con)

}

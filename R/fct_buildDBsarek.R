#' buildDB_sarek
#' @description
#' A short description...
#' @param name description
#' @importFrom VariantAnnotation readVcf scanVcfHeader geno header info
#' @importFrom tibble tibble
#' @importFrom magrittr %<>% 
#' @importFrom dplyr %>% select filter mutate inner_join group_by case_when bind_cols mutate_all mutate_if arrange as_tibble
#' @importFrom stringr str_replace_all str_detect str_c str_split
#' @importFrom tidyr separate unnest gather spread replace_na pivot_wider
#' @importFrom rlang enquo
#' @importFrom reshape2 melt
#' @importFrom IRanges stack
#' @importFrom promises promise_all future_promise %...>%
#' @return The return value, if any, from executing the function.
#'
#' @noRd
#' @export
buildDB_sarek <- function(prefix = NULL, vcf_name = NULL, db_path = NULL) {
  
  ############### uncomment for dev tests #######################################
  #db_path <- system.file("extdata","testdata", package = "SomaVarDB")
  #prefix <- "test"
  #vcf_name <- system.file("extdata","testdata/vcf_name", package = "SomaVarDB")
  #library(VariantAnnotation)
  #library(tibble)  
  #library(magrittr) 
  #library(dplyr)
  #library(stringr)
  #library(tidyr)
  #library(rlang)
  #library(promises)
  #library(DBI)
  #library(RSQLite)
  ##############################################################################
  
  process_sample <- FALSE
  db_name <- file.path(db_path,paste0(prefix, ".db"))
  #### A function to convert genotypes to dosage
  gt2snp <- function(x){
    case_when(
      str_detect(x, "0\\/\\.|\\.\\/0|0\\/0|0\\|0|^0$")  ~ 0,
      str_detect(x, c("0\\/1|1\\/0|0\\|1|1\\|0|\\.\\/1|1\\/\\.|\\.\\|1|1\\|\\.|^1$"))  ~ 1,
      str_detect(x,  c("1\\/1|1\\|1")) ~ 2,
      TRUE ~ as.numeric(NA))
  }
  #### read VCF information ####
  message('######\nPARSING VCF HEADER\n######')
  vcf_header <- scanVcfHeader(vcf_name)
  .vcf <- readVcf(vcf_name)
  samples_frame <- data.frame(sample = vcf_header@samples)
  samples_frame$constant_group <- "ALL_DB"
  
  #### Is sample already in the db ####
  con <- dbConnect(SQLite(), db_name)
  if(DBI::dbExistsTable(conn = con, name = "samples")){
    samples <- DBI::dbReadTable(conn = con,name = "samples")
    if (!(vcf_header@samples %in% samples$sample)){
      samples_frame <- bind_rows(samples,samples_frame)
      DBI::dbWriteTable(conn = con, name = "samples", value= samples_frame, overwrite = TRUE)
      process_sample <- TRUE}
  } else {
    DBI::dbWriteTable(conn = con, name = 'samples',value = samples_frame)
    process_sample <- TRUE
  }
  DBI::dbDisconnect(con)
  
  if(process_sample){
    #### Check to see if CSQ exists ####
    if(class(.vcf@info$CSQ)!="NULL"){
      csq_exists = TRUE
      csq_cols <- vcf_header@header$INFO$Description[rownames(vcf_header@header$INFO) == 'CSQ'] %>%
        tolower() %>% str_replace_all(" |:|\\.", "_")
    } else  {
      csq_exists = FALSE
      message("######\nNo VEP annotations detected.\n######")
    }
    ## These are the impacts flagged as 'exonic' from VEP
    exonic_impacts = c("stop_gained","exon_variant",
                       "stop_lost","frameshift_variant","initiator_codon_variant",
                       "inframe_deletion","inframe_insertion",
                       "missense_variant","protein_altering_variant","incomplete_terminal_codon_variant",
                       "stop_retained_variant","5_prime_UTR_premature_start_codon_variant",
                       "synonymous_variant","coding_sequence_variant",
                       "5_prime_UTR_variant", "3_prime_UTR_variant", "transcript_ablation",
                       "transcript_amplification","feature_elongation","feature_truncation")
    
    #### Build Tables ####
    message("######\nSTARTING TO BUILD DATABASE\n######")
    vcf_info <- .vcf@info %>% as_tibble()
    vcf_rowRanges  <- .vcf@rowRanges %>% as_tibble()
    vcf_fixed  <- .vcf@fixed %>% as_tibble()
    vcf_geno <- geno(.vcf)
    names <- names(.vcf)
    rm(.vcf)
    start <- as.numeric(gsub("_.*","",gsub(".*:","",names)))
    mut <- gsub(".*_", "", names)
    var_ind <- paste0((start - 1),"_",mut)
    
    csq_promise <- future_promise({
      if(csq_exists){
        csq.vcf <- vcf_info %>%
          dplyr::select(CSQ) %>%
          tibble::add_column(variant_id = var_ind) %>%
          unnest(CSQ) %>%
          separate(CSQ, sep = "\\|", into = str_split(csq_cols, pattern = "\\|")[[1]]) %>%
          #separate_rows(consequence,  sep = "&") %>%
          mutate(is_lof = impact == "HIGH" & biotype == "protein_coding",
                 is_splicing = str_detect(consequence, "splice"),
                 is_exonic = biotype %in% exonic_impacts,
                 is_intronic = intron != "",
                 VKB = "Unknown", commentary = "") %>%
          dplyr::select(-ends_with("_af"), -any_of(c('clin_sig', 'pheno','somatic', 'pubmed', 'consequence_annotations_from_ensembl_vep__format__allele')))
        
        csq.vcf$hgvsp_vep <- gsub("^.*:","",csq.vcf$hgvsp)
        
        amino_acid_mapping <- list(
          Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
          Glu = "E", Gln = "Q", Gly = "G", His = "H", Ile = "I",
          Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
          Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V"
        )
        
        # Function to convert three-letter amino acid code to single-letter code
        convert_to_single_letter <- function(protein_change) {
          for (three_letter in names(amino_acid_mapping)) {
            single_letter <- amino_acid_mapping[[three_letter]]
            protein_change <- gsub(three_letter, single_letter, protein_change)
            print(single_letter)
            print(protein_change)
          }
          return(protein_change)
        }
        
        csq.vcf$hgvsp <- sapply(csq.vcf$hgvsp_vep, function(line) {
          single_letter_change <- convert_to_single_letter(line)
          return(single_letter_change)
        })
        
        clinvar.vcf <- vcf_info %>%
          dplyr::select(any_of(c("clinvar_sig", "clinvar_disease_name"))) %>%
          as_tibble() %>%
          bind_cols(tibble(variant_id = var_ind))
        
        csq.vcf <- csq.vcf %>%
          {if('clinvar_sig' %in% colnames(clinvar.vcf)){
            left_join(., clinvar.vcf %>%
                        dplyr::select(variant_id, clinvar_sig) %>%
                        unnest(clinvar_sig),
                      by = "variant_id")
          } else {.}} %>%
          {if('clinvar_disease_name' %in% colnames(clinvar.vcf)){
            left_join(., clinvar.vcf %>%
                        dplyr::select(variant_id, clinvar_disease_name) %>%
                        unnest(clinvar_disease_name),
                      by = "variant_id")
          } else {.}} %>%
          dplyr::select(variant_id, everything())
        return(csq.vcf)
      }
    }) # end of csq future promises
    
    info_promise <- future_promise({
      if(csq_exists){
        info.vcf <- tibble(variant_id = var_ind) %>%
          bind_cols(vcf_rowRanges  %>% mutate(chr = as.character(seqnames)) %>% dplyr::select(chr, start, end)) %>%
          bind_cols(vcf_fixed) %>%
          bind_cols(vcf_info %>%  dplyr::select(-any_of(c('CSQ', 'clinvar_sig', 'clinvar_disease_name'))))
      } else {
        print("CSQ DOES NOT EXISTS")
        info.vcf <- tibble(variant_id = var_ind) %>%
          bind_cols(vcf_rowRanges  %>%
                      dplyr::select(seqnames, start, end) %>%
                      rename('chr' = seqnames)) %>%
          bind_cols(vcf_fixed) %>%
          bind_cols(vcf_info)
      }
      names(info.vcf) %<>% tolower()
      info.vcf <- info.vcf[,!duplicated(colnames(info.vcf))]
      info.vcf <- info.vcf %>%
        mutate(alt = sapply(alt, as.character)) %>%
        mutate_if(function(x){class(x)=="AsIs"}, as.character) %>%
        #dplyr::select(all_of(c("variant_id","chr","start","end","ref","alt","qual","filter","ab","abp","ac","af","an","ao","cigar","min_dp","mqm",
        #                       "mqmr","ns","numalt","odds","ro","rpp","rppr","rpr","run","saf","sap","sar","type"))) %>% # variant info to keep for the moment. Has to be common accross all vcf sources
        mutate_all(as.character)
    }) # end of info promise
    
    ## fix a few weird formatting things on the columns
    ## then add filepaths for variants
    geno_col_names <- vcf_geno %>% names()
    geno_promise <- future_promise({
      if (length(geno_col_names) > 0){
        .geno_col <- geno_col_names[1]
        geno_col <- enquo(.geno_col)
        geno.vcf <- tibble(group = var_ind, variant_id = var_ind) %>%
          bind_cols(vcf_geno[[.geno_col]] %>% as_tibble()) %>%
          gather(sample, !!.geno_col, -variant_id, -group)
        
        for(.geno_col in geno_col_names[-1]){
          message(.geno_col)
          geno_col <- enquo(.geno_col)
          if ("matrix" %in% class(vcf_geno[[.geno_col]])){
            geno.vcf <- geno.vcf %>%
              bind_cols(
                vcf_geno[[.geno_col]] %>%
                  as_tibble() %>% gather(sample, !!.geno_col) %>%
                  dplyr::select(-sample))
          } else {
            tmp <- vcf_geno[[.geno_col]]
            rownames(tmp) <- var_ind
            tmp <- tmp %>%
              reshape2::melt() %>%
              dplyr::select(variant_id = Var1, sample = Var2, Var3, value) %>%
              mutate(Var3 = str_c(.geno_col, "_", Var3)) %>%
              spread(Var3, value)
            
            geno.vcf <- geno.vcf %>%
              left_join(tmp, by = c("variant_id", "sample"))
          }
        }
        names(geno.vcf) %<>% tolower()
        
        if('GT' %in% geno_col_names){
          geno.vcf <- geno.vcf %>%
            mutate(gt_raw = gt, gt = gt2snp(gt_raw))
        }
        geno.vcf <- geno.vcf %>%
          dplyr::select(-group) %>%
          arrange(variant_id, sample) %>%
          ##### SEQONE VERSION #########
        #mutate(ao = as.numeric(ao),
        #       af = round(ao/as.numeric(dp), digits = 2), # because empty on seqone results
        #       qa = 10*log10(gq)) %>%
        ## gq = probability that the call is incorrect = genotype quality
        ##  qa = same as phred score
        #########################         
        #      filter(qa >= 10)
        ##### SAREK VERSION #########
        separate(ad, c("delete","ro","ao")) %>%
          mutate(#ao = as.numeric(ao),
            af = round(as.numeric(ao)/as.numeric(dp), digits = 2), # because empty on seqone results
            qa = 10*log10(as.numeric(gq))) %>%
          ## gq = probability that the call is incorrect = genotype quality
          ##  qa = same as phred score
          #########################         
        filter(qa >= 10)
        geno.vcf <- geno.vcf[,!(colnames(geno.vcf) %in% c("gq","pl","delete"))]
      }
    }) # end of future promises geno
    
    promise_all <- promise_all(geno_promise, csq_promise, info_promise) %...>% {
      
      geno.vcf <- environment(geno_promise[["then"]])[["private"]][["value"]]
      info.vcf <- environment(info_promise[["then"]])[["private"]][["value"]]
      csq.vcf <- environment(csq_promise[["then"]])[["private"]][["value"]]
      
      message("######\n Inserting vcf metadata \n#####")
      #header <- vcf_header@header[c("AltaiNeandertal","Denisova","GeneSplicer","VEP","aapos","codonpos","commandline","fileDate","fileformat","phasing")]
      #### QUELLES METADATA A DEFINIR DANS LE HEADER AVANT LA BASE #### AU DESSUS SONT CELLES VIA VCF SEQONE
      
      header <- data.frame("AltaiNeandertal" = "sarek",
                           "Denisova" = "sarek",
                           "GeneSplicer" = "sarek",
                           "VEP" = "sarek",
                           "aapos" = "sarek",
                           "codonpos" = "sarek",
                           "commandline" = "sarek",
                           "fileDate" = "sarek",
                           "fileformat" = "sarek",
                           "phasing" = "sarek" )
      
      vcf_metadata <- as.data.frame(IRanges::stack(header, index.var = "names"))
      vcf_metadata$sample <- vcf_header@samples
      vcf_metadata <- vcf_metadata %>% pivot_wider(names_from = "ind", values_from = "values")
      con <- dbConnect(SQLite(), db_name)
      DBI::dbWriteTable(conn = con, name = "vcf_metadata", value = vcf_metadata, append = TRUE)
      
      message("######\nStarting inserting variants\n#####")
      if(DBI::dbExistsTable(conn = con,name="variant_impact")){
        DBI::dbWriteTable(conn = con,name="variant_impact_tmp", value =  csq.vcf, overwrite = TRUE)
        dbSendQuery(con,"CREATE INDEX idx_variant_impact_tmp ON variant_impact_tmp (variant_id);")
        
        # buildSeqonemethod :
        # apn_sql <- "INSERT INTO variant_impact SELECT * FROM variant_impact_tmp s WHERE NOT EXISTS (SELECT 1 FROM variant_impact t WHERE t.[variant_id] = s.[variant_id]);"
        
        
        # Try to implement mixed sources compatibiliy on variant_impact table, keeping only commond columns
        # not_in_tmp <- dbGetQuery(con,"SELECT name FROM PRAGMA_TABLE_INFO('variant_impact') WHERE name NOT IN (SELECT name FROM PRAGMA_TABLE_INFO('variant_impact_tmp'));")
        # impact_colnames <- dbGetQuery(con,"PRAGMA name table_info(variant_impact);")
        # paste0("INSERT INTO variant_impact (", paste0(impact_colnames$name,collapse = ","),") SELECT ",paste0(common$common_colnames,collapse = ","), "")
        common <- dbGetQuery(con,"SELECT variant_impact.name as common_colnames FROM pragma_table_info('variant_impact') AS variant_impact JOIN pragma_table_info('variant_impact_tmp') AS variant_impact_tmp ON variant_impact.name  = variant_impact_tmp.name;")
        dbSendQuery(con, paste0("INSERT INTO variant_impact (", paste0(common$common_colnames,collapse = ","),") SELECT ",paste0(common$common_colnames,collapse = ","), " FROM variant_impact_tmp s WHERE NOT EXISTS (SELECT 1 FROM variant_info t WHERE t.[variant_id] = s.[variant_id]);"))
        try({dbSendQuery(con,"DROP INDEX idx_variant_impact;")})
        dbSendQuery(con,"CREATE INDEX idx_variant_impact ON variant_impact (variant_id);")
        
      } else {
        DBI::dbWriteTable(
          conn = con, name = "variant_impact",value = csq.vcf, append = TRUE)
        dbSendQuery(con,"CREATE INDEX idx_variant_impact ON variant_impact (variant_id);")
      }
      
      if(DBI::dbExistsTable(conn = con, name="variant_info")){
        DBI::dbWriteTable( conn = con, name = "variant_info_tmp", value = info.vcf, overwrite = TRUE)
        dbSendQuery(con,"CREATE INDEX idx_variant_info_tmp ON variant_info_tmp (variant_id);")
        
        # buildSeqonemethod :
        #apn_sql <- "INSERT INTO variant_info SELECT * FROM variant_info_tmp s WHERE NOT EXISTS (SELECT 1 FROM variant_info t WHERE t.[variant_id] = s.[variant_id]);"
        #dbSendQuery(con, apn_sql)
        
        # Try to implement mixed sources compatibiliy on variant_info table, keeping only commond columns
        common <- dbGetQuery(con,"SELECT variant_info.name as common_colnames FROM pragma_table_info('variant_info') AS variant_info JOIN pragma_table_info('variant_info_tmp') AS variant_info_tmp ON variant_info.name  = variant_info_tmp.name;")
        dbSendQuery(con, paste0("INSERT INTO variant_info (", paste0(common$common_colnames,collapse = ","),") SELECT ",paste0(common$common_colnames,collapse = ","), " FROM variant_info_tmp s WHERE NOT EXISTS (SELECT 1 FROM variant_info t WHERE t.[variant_id] = s.[variant_id]);"))
        try({dbSendQuery(con,"DROP INDEX idx_variant_info;")})
        dbSendQuery(con,"CREATE INDEX idx_variant_info ON variant_info (variant_id);")
        
      } else {
        DBI::dbWriteTable( conn = con, name = "variant_info", value = info.vcf , append = TRUE)
        dbSendQuery(con,"CREATE INDEX idx_variant_info ON variant_info (variant_id);")
      }
      
      # buildSeqonemethod :
      # DBI::dbWriteTable(conn = con, name = "variant_geno", value = geno.vcf, append = TRUE)
      # try({dbSendQuery(con,"DROP INDEX idx_variant_geno;")})
      # dbSendQuery(con,"CREATE INDEX idx_variant_geno ON variant_geno (variant_id);")
      
      
      # different sources compatibilities 
      if(DBI::dbExistsTable(conn = con, name="variant_geno")){
        DBI::dbWriteTable(conn = con, name = "variant_geno_tmp", value = geno.vcf, overwrite = TRUE)
        common <- dbGetQuery(con,"SELECT variant_geno.name as common_colnames FROM pragma_table_info('variant_geno') AS variant_geno JOIN pragma_table_info('variant_geno_tmp') AS variant_geno_tmp ON variant_geno.name  = variant_geno_tmp.name;")
        dbSendQuery(con, paste0("INSERT INTO variant_geno (", paste0(common$common_colnames,collapse = ","),") SELECT ",paste0(common$common_colnames,collapse = ","), " FROM variant_geno_tmp;"))
        try({dbSendQuery(con,"DROP INDEX idx_variant_geno;")})
        dbSendQuery(con,"CREATE INDEX idx_variant_geno ON variant_geno (variant_id);")
      } else {
        DBI::dbWriteTable( conn = con, name = "variant_geno", value = geno.vcf , append = TRUE)
        dbSendQuery(con,"CREATE INDEX idx_geno_info ON variant_geno (variant_id);")
      }
      
      DBI::dbDisconnect(con)
      closeAllConnections()
      message("######\nDone inserting variants\n#####")
    } # end or promise concatenation
    
    block_until_settled <- function(p) {
      promise_resolved <- FALSE
      p$finally(function(value) {
        promise_resolved <<- TRUE
      })
      
      while(!promise_resolved) {
        later::run_now(1)
      }
    }
    
    print(promise_all)
    
    block_until_settled(promise_all)
    
  } else {cat("sample already present in the database, skipping import")}
  gc(verbose = TRUE)
}

#' vcftodb
#'
#' @description A fct function
#' @importFrom jsonlite fromJSON
#' @importFrom httr POST GET
#' @importFrom RSQLite SQLite
#' @importFrom DBI dbGetQuery dbExistsTable dbWriteTable dbDisconnect dbConnect
#' @importFrom crew crew_controller_local
#' @importFrom dplyr %>% select mutate rename slice
#'
#' @noRd
#' @export
addMDtodb <- function(API_key = "",
                      prefix = NULL,
                      db_path = NULL, workers = 2) {
  
  if(API_key %in% c(NULL,"")){
   print("A Mobidetails API key must be provided to the API_key argument...")
   print("To obtain an API key create a MobiDetails account : https://mobidetails.iurc.montp.inserm.fr/MD ")
   if(readline("Press q to exit the code: ") == 'q'){
     .Internal(.invokeRestart(list(NULL,NULL), NULL))
   }
  return(1)
  }
  
  options("HTTPUserAgent" = "CHUGA")
  db_name <- file.path(db_path,paste0(prefix, ".db"))
  con <- dbConnect(SQLite(), db_name)
  controller <- crew_controller_local(name = "Mobideltails in parallel ", workers = workers,  seconds_idle = 100)
  
  variants_info <- data.frame(notempty = "yes")
  while(nrow(variants_info) > 0){
    already_annoted <- if(DBI::dbExistsTable(con,"variant_MD")){dbGetQuery(con,"SELECT variant_id FROM variant_MD;")}else{data.frame(variant_id = "variant0")}
    variants_info <- dbGetQuery(con, paste0("SELECT * from variant_info WHERE variant_id NOT IN ('",
                                            paste0(already_annoted$variant_id,collapse="' , '"),
                                            "');")) %>% dplyr::slice(1:100)

    if(nrow(variants_info) == 0){print("end of annotation process");try({controller$terminate()});break}
    else {
      controller$start()
      variants_MD <- controller$map(
        command = tryCatch(
          expr = {
            row <- variants_info[nrow,]
            variant <- row[,"variant_id"]
            print(paste0("processing variant : ", variant))
            variant_query <- paste0(row[,"chr"],"-",row[,"start"] ,"-",row[,"ref"],"-",row[,"alt"])
            url <- paste0("https://mobidetails.iurc.montp.inserm.fr/MD/api/variant/create_vcf_str?genome_version=hg19&vcf_str=",variant_query,"&caller=cli&api_key=",API_key)
            raw_data <- httr::POST(url,httr::verbose(),encode="json",
                                   httr::user_agent(getOption("HTTPUserAgent")))
            json <- jsonlite::fromJSON(rawToChar(raw_data$content))
            mdurl <- paste0('<a href="',json[[1]]$url,'"target="_blank" class="btn btn-primary"',">",variant,"</a>")
            print(mdurl)

            url <- paste0("https://mobidetails.iurc.montp.inserm.fr/MD/api/variant/",json[[1]]$mobidetails_id,"/cli/,")
            raw_data <- httr::GET(url,httr::verbose(),encode="json",
                                httr::user_agent(getOption("HTTPUserAgent")))
            json <- jsonlite::fromJSON(rawToChar(raw_data$content))
            frame <- as.data.frame(t(as.data.frame(lapply(json, function(x){
              x[sapply(x,is.null)] <- NA
              unlist(unlist(x))
            }) %>% unlist())))
            rownames(frame) <- variant
            frame <- frame %>% dplyr::select(c(
              "gene.isTumorSuppressor",
              "gene.isOncogene",
              #"gene.canonical",
              "frequenciesDatabases.dbSNPrsid",
              "frequenciesDatabases.gnomADv3",
              "missensePredictions.siftPred"  ,
              "missensePredictions.siftScore" ,
              "missensePredictions.polyphen2HdivPred" ,
              "missensePredictions.polyphen2HdivScore" ,
              "missensePredictions.polyphen2HvarPred" ,
              "missensePredictions.polyphen2HvarScore" ,
              "frequenciesDatabases.clinvarClinsig" ,
              "frequenciesDatabases.clinvarClinsigConf"))
            frame["variant_id"] <- variant
            frame["mdurl"] <- mdurl
            frame["isMDetailed"]  <- "Yes"
            frame <- frame %>% dplyr::select("variant_id",everything())
            return(frame)
          }, 
          error = function(e){
            print(paste("MDetails error for : ", variant))
            print(e)
            frame <- data.frame("variant_id" = variant,
                                "mdurl" = variant,
                                "isMDetailed" = "Absent on MobiDetails",
                                "gene.isTumorSuppressor" = "Absent on MobiDetails",
                                "gene.isOncogene" = "Absent on MobiDetails",
                                #"gene.canonical" = "Absent on MobiDetails",
                                "frequenciesDatabases.dbSNPrsid" = "Absent on MobiDetails",
                                "frequenciesDatabases.gnomADv3" = "Absent on MobiDetails",
                                "missensePredictions.siftPred"  = "Absent on MobiDetails",
                                "missensePredictions.siftScore" = "Absent on MobiDetails",
                                "missensePredictions.polyphen2HdivPred" = "Absent on MobiDetails",
                                "missensePredictions.polyphen2HdivScore" = "Absent on MobiDetails",
                                "missensePredictions.polyphen2HvarPred" = "Absent on MobiDetails",
                                "missensePredictions.polyphen2HvarScore" = "Absent on MobiDetails",
                                "frequenciesDatabases.clinvarClinsig" = "Absent on MobiDetails",
                                "frequenciesDatabases.clinvarClinsigConf" = "Absent on MobiDetails")                             
            rownames(frame) <- variant
            return(frame)
      }),
        iterate = list(nrow = seq(1:nrow(variants_info))),
        globals = list(variants_info = variants_info, API_key = API_key),
        packages = c("dplyr")
      )
      controller$terminate()

      variants_MD <- do.call(rbind,variants_MD$result)
      #return(variants_MD$result)
      variants_MD <- variants_MD %>% dplyr::rename("TumorSuppressor" = "gene.isTumorSuppressor",
                                                  "Oncogene" = "gene.isOncogene", 
                                                  #"canonical" = "gene.canonical",
                                                  "dbSNP" = "frequenciesDatabases.dbSNPrsid",
                                                  "gnomADv3" = "frequenciesDatabases.gnomADv3",
                                                  "siftPred" = "missensePredictions.siftPred",
                                                  "siftScore" = "missensePredictions.siftScore",
                                                  "polyphen2HdivPred" = "missensePredictions.polyphen2HdivPred",
                                                  "polyphen2HdivScore" = "missensePredictions.polyphen2HdivScore",
                                                  "polyphen2HvarPred" = "missensePredictions.polyphen2HvarPred",
                                                  "polyphen2HvarScore" = "missensePredictions.polyphen2HvarScore",
                                                  "clinvarClinsig" = "frequenciesDatabases.clinvarClinsig",
                                                  "clinvarClinsigConf" = "frequenciesDatabases.clinvarClinsigConf")

      if(DBI::dbExistsTable(con,"variant_MD")){DBI::dbWriteTable(conn = con, name = "variant_MD", value = variants_MD, append = TRUE)}
      else{DBI::dbWriteTable(conn = con, name = "variant_MD", value = variants_MD, overwrite = TRUE)}
    }
  }
  DBI::dbDisconnect(con)
}

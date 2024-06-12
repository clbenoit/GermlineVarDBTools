#' Download reference genome for local JBrowseR server
#'
#' @param ... arguments to pass to the function
#' params description
#'
#' @export
download_genome_reference <- function(browser_server_path, genome_version) {
  print("entering jbrowser local set up function")
  if(genome_version == "hg19"){
    if(!file.exists(file.path(browser_server_path, "human_g1k_v37.fasta.gz"))){
       system(command = paste0("cd ",browser_server_path, 
                               "; wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
                               )
              )
    } else {print("human_g1k_v37.fasta.gz already present, skipping download")}
  } else if (genome_version == "hg38"){
    if(!file.exists(file.path(browser_server_path, "GRCh38_full_analysis_set_plus_decoy_hla.fa"))){
      system(command = paste0("cd ",
                              browser_server_path, 
                              "; wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
                              )
             )
    } else {print("GRCh38_full_analysis_set_plus_decoy_hla.fa already present, skipping download")}
  }
}
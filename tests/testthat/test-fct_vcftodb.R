testthat::test_that("vcftodb works", {
  
  db_path <- system.file("extdata","testdata", package = "GermlineVarDBTools")
  prefix <- "test"
  vcf_list <- list.files(system.file("extdata","testdata", package = "GermlineVarDBTools"), pattern = ".gz$", full.names = TRUE)

  future::plan("multicore")
  options(future.globals.maxSize = 10000*1024^2)
  options(mc.cores = 3)
  for (vcf in vcf_list){
    print(vcf)
    GermlineVarDBTools::buildDB_seqone(prefix = prefix, vcf_name = vcf, db_path = db_path)
  }

  GermlineVarDBTools::addMDtodb(db_path = db_path, prefix = prefix, API_key = "l6ln_tsCbdZk8ICVpy6z86vhhb2rFUkFFTfV3wq7L-4", workers = 2)
  
  GermlineVarDBTools::update_db_metadata(db_path = db_path, prefix = prefix)
  
  GermlineVarDBTools::addSampleAttributes(db_path = db_path, prefix = prefix, matches_file = system.file("extdata","testdata/samples_run_match.tsv", package = "GermlineVarDBTools"))

  GermlineVarDBTools::compute_frequency(db_path = db_path, prefix = prefix ,attribute = "constant_group")
  
  GermlineVarDBTools::CNVtodb_seqone(db_path = db_path, prefix = prefix, cnvfile_path = system.file("extdata","testdata/cnv.tsv", package = "GermlineVarDBTools"))

  GermlineVarDBTools::QCtodb_seqone(db_path = db_path, prefix = prefix, qcfile_path = system.file("extdata","testdata/QCs.csv", package = "GermlineVarDBTools"))
  
  GermlineVarDBTools::initcomment(db_path = db_path, prefix = prefix)
})

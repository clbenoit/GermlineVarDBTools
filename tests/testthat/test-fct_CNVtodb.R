testthat::test_that("CNVtodb works", {
  cnvfile_path = system.file("extdata","testdata/cnv.tsv", package = "GermlineVarDB")
  prefix = "test"
  db_path = system.file("extdata","testdata", package = "GermlineVarDB")
  CNVtodb_seqone(cnvfile_path = cnvfile_path,db_path =db_path,
                   prefix = prefix)
})

#### CHECKER SI LES DATA SONT DUPLIQUUES PAR SAMPLE/GENE ET PK ?

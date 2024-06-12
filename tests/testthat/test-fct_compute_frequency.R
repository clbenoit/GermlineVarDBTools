testthat::test_that("frequency computation works", {
  db_path = system.file("extdata","testdata", package = "GermlineVarDB")
  prefix = "test"
  GermlineVarDB::compute_frequency(db_path = db_path, prefix = prefix,attribute = "constant_group")
})
  

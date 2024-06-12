test_that("database metadata updating works", {
  update_db_metadata(
          db_path = system.file("extdata","testdata", package = "GermlineVarDB"),
          prefix = "test")

})

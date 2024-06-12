testthat::test_that("adding Mobidetails to database works", {
  addMDtodb(API_key = "l6ln_tsCbdZk8ICVpy6z86vhhb2rFUkFFTfV3wq7L-4",
            #db_path = dirname(system.file("extdata","dbs/test.db", package = "GermlineVarDB")),
            db_path = system.file("extdata","testdata", package = "GermlineVarDB"),
            prefix = "test")
})

test_that("bac.df is correctly loaded", {
  data(bac.df)
  expect_true(is.data.frame(bac.df))
  expect_named(bac.df, c("POS", "QUAL", "ISO", "GQ", "REF", "ALT", "COV", "COV.REF",
                         "COV.ALT1", "COV.ALT2", "COV.ALT3", "REL.REF", "REL.ALT1", "REL.ALT2",
                         "REL.ALT3"))
})

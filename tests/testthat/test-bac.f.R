test_that("bac.f is correctly loaded", {
  data(bac.f)
  expect_true(is.data.frame(bac.f))
  expect_named(bac.f, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                        "GT", "PL", "DP", "DPR", "GQ", "ISO", "COV", "COV.REF", "COV.ALT1",
                        "REL.REF", "REL.ALT1"))
})

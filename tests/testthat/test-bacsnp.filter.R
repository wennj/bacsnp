test_that("bacsnp.filter filters correctly", {
  data(bac.df)
  # Minimum test
  filtered <- bacsnp.filter(bac.df, min.abs.cov = 10)
  expect_true(is.data.frame(filtered))
  expect_true(nrow(filtered) <= nrow(bac.df))  # No additional rows added

  # Regular test
  filtered2 <- bacsnp.filter(bac.df, min.abs.cov = 100, min.abs.alt = 10, min.rel.alt = 0.05)
  expect_true(nrow(filtered2) <= nrow(bac.df))
})

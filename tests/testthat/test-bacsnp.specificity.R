test_that("bacsnp.specificity calculates correctly", {
  data(bac.f)
  result <- bacsnp.specificity(bac.f, isolates = c("iIso1", "iIso2"), which.rel = "REL.ALT1")
  expect_true(is.list(result))
  expect_true(length(result) == 2)  # Sollte 2 Elemente enthalten
})

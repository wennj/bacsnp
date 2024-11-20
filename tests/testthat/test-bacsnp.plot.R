test_that("bacsnp.plot generates ggplot object", {
  data(bac.f)
  plot <- bacsnp.plot(bac.f, col = "ISO", genome.length = 123193, mark.repeats = FALSE)
  expect_s3_class(plot, "ggplot")
})

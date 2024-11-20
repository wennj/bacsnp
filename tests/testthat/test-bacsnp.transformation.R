test_that("bacsnp.transformation converts VCF to dataframe", {
  # real vcf file as input
  vcf_example <- vcfR::read.vcfR(system.file("extdata", "bac.vcf", package = "bacsnp"))
  result <- bacsnp.transformation(vcf_example)
  expect_true(is.data.frame(result))
  expect_gt(ncol(result), 1)  # Should have more than one column
})

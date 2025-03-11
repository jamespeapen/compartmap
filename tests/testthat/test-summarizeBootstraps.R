# summarizeBoostraps
test_that("summarizeBoostraps", {
})

# .isCompartmentOpen
test_that(".isCompartmentOpen", {
  .isCompartmentOpen <- compartmap:::.isCompartmentOpen
  expect_true(
    .isCompartmentOpen(is.atac_or_rna = TRUE, eigen = 1),
    .isCompartmentOpen(is.atac_or_rna = FALSE, eigen = 0),
    .isCompartmentOpen(is.atac_or_rna = FALSE, eigen = -1)
  )
  expect_false(
    .isCompartmentOpen(is.atac_or_rna = TRUE, eigen = 0),
    .isCompartmentOpen(is.atac_or_rna = TRUE, eigen = -1),
    .isCompartmentOpen(is.atac_or_rna = FALSE, eigen = 1)
  )
})

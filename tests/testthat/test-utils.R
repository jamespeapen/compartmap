library(SummarizedExperiment)
library(RaggedExperiment)
library(GenomicRanges)

se <- SummarizedExperiment()
re <- RaggedExperiment()
df <- data.frame()
chrs <- c("chr2", "chr2", "chr1", "chr3")
gr <- GRanges(
  Rle(chrs, c(1, 3, 2, 4)),
  IRanges(1:10, width=10:1)
)

test_that("checkAssayType", {
  expect_true(checkAssayType(se))
  expect_false(checkAssayType(re))
  expect_false(checkAssayType(gr))
  expect_false(checkAssayType(df))
})


gr.pc <- gr
pc.list <- list('pc' =  c(1,        -1,       1,         -1,       -2,       2,       -5,        5,      -10,     10))
pc.state.cutoff.zero <- c("open",   "closed", "open",   "closed", "closed", "open",   "closed", "open", "closed", "open")
pc.state.cutoff.two <-  c("closed", "closed", "closed", "closed", "closed", "closed", "closed", "open", "closed", "open")
pc.state.array <-  lapply(pc.state.cutoff.zero, function(i) ifelse(i == "open", "closed", "open")) |> unlist()
mcols(gr.pc) <- pc.list

test_that("extractOpenClosed", {
  expect_error(extractOpenClosed(se))
  expect_error(extractOpenClosed(re))
  expect_error(extractOpenClosed(df))
  expect_error(extractOpenClosed(gr))
  expect_no_error(extractOpenClosed(gr.pc))
  expect_equal(
    extractOpenClosed(gr.pc, assay = "rna"),
    pc.state.cutoff.zero
  )
  expect_equal(
    extractOpenClosed(gr.pc, assay = "rna", cutoff = 2),
    pc.state.cutoff.two
  )
  expect_equal(
    extractOpenClosed(gr.pc, assay = "atac"),
    pc.state.cutoff.zero
  )
  expect_equal(
    extractOpenClosed(gr.pc, assay = "array"),
    pc.state.array
  )
})

test_that("getChrs", {
  expect_equal(
    getChrs(gr),
    unique(chrs)
  )
})


nrows <- 200
ncols <- 6
counts1 <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
counts2 <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
counts3 <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
colData <- DataFrame(
  Treatment = rep(c("a", "b"), 3),
  row.names = LETTERS[1:6]
)
se <- SummarizedExperiment(
  assays = SimpleList(one = counts1, two = counts2, three = counts3),
  colData = colData
)

test_that("getAssayNames", {
  expect_equal(
    getAssayNames(se),
    c("one", "two", "three")
  )
})

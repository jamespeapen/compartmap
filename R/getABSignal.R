#' Calculate Pearson correlations of smoothed eigenvectors
#'
#' This function is used to generate a list x to be passed to getABSignal
#'
#' @param x      A list object from getCorMatrix
#' @param squeeze    Whether squeezing was used (implies Fisher's Z transformation)
#' @param assay What kind of assay are we working on ("array", "atac", "array")
#' @param genome The genome to use for gene-density-based sign correction
#' @param smooth Whether to smooth the singular vector
#'
#' @return    A list x to pass to getABSignal
#'
#' @import    SummarizedExperiment
#'
#' @export
#'
#' @examples
#'
#' library(SummarizedExperiment)
#' library(BiocSingular)
#'
#' #Generate random genomic intervals of 1-1000 bp on chr1-22
#' #Modified from https://www.biostars.org/p/225520/
#' random_genomic_int <- data.frame(chr = rep("chr14", 100))
#' random_genomic_int$start <- apply(random_genomic_int, 1, function(x) {
#'   round(runif(1, 0, getSeqLengths(getGenome("hg19"), chr = x)[[1]]), 0)
#' })

#' random_genomic_int$end <- random_genomic_int$start + runif(1, 1, 1000)
#' random_genomic_int$strand <- "*"
#'
#' #Generate random counts
#' counts <- rnbinom(1000, 1.2, 0.4)
#'
#' #Build random counts for 10 samples
#' count.mat <- matrix(sample(counts, nrow(random_genomic_int) * 10, replace = FALSE), ncol = 10)
#' colnames(count.mat) <- paste0("sample_", seq(1:10))
#'
#' #Bin counts
#' bin.counts <- getBinMatrix(
#'   count.mat,
#'   makeGRangesFromDataFrame(random_genomic_int),
#'   chr = "chr14",
#'   genome = "hg19"
#' )
#'
#' #Calculate correlations
#' bin.cor.counts <- getCorMatrix(bin.counts)
#'
#' #Get A/B signal
#' absignal <- getABSignal(bin.cor.counts)

getABSignal <- function(
  x,
  squeeze = FALSE,
  assay = c("rna", "atac", "array"),
  genome = c("hg19", "hg38", "mm9", "mm10"),
  smooth = TRUE
) {
  assay <- match.arg(assay)
  gen <- match.arg(genome)
  gr <- x$gr

  flog.debug("Calculating eigenvectors.")
  cscore <- getSVD(x$binmat.cor, sing.vec = "right") |> as.vector()
  if (squeeze) {
    cscore <- ifisherZ(cscore)
  }

  if (smooth) {
    flog.debug("Smoothing eigenvector.")
    gr$cscore <- meanSmoother(cscore)
    flog.debug("Done smoothing.")
  } else {
    gr$cscore <- cscore
  }

  if (flipSign(gr, genome, assay)) {
    gr$cscore <- -gr$cscore
  }
  Seqinfo::genome(gr) <- gen
  return(gr)
}

flipSign <- function(gr, genome, assay) {
  tx.gr <- getGenome(genome, "tx")
  gene_count <- countOverlaps(gr, tx.gr)
  open <- gr$cscore > 0
  flip <- sum(gene_count[open]) < sum(gene_count[!open])
  flip
}

# Check if a compartment is open based on assay type and eigenvalue
#
# For ATAC/RNA:
# eigen < cutoff - closed
# eigen > cutoff - open
#
# For methylation the logic is flipped:
# eigen < cutoff - open
# eigen > cutoff - closed
.isCompartmentOpen <- function(is.atac_or_rna, eigen, cutoff) {
  (is.atac_or_rna & eigen > cutoff) | (!is.atac_or_rna & eigen < cutoff)
}

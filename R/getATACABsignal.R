#' @title Estimate A/B compartments from ATAC-seq data
#'
#' @description
#' \code{getATACABsignal} returns estimated A/B compartments from ATAC-seq data.
#'
#' @param obj Input SummarizedExperiment object
#' @param res Compartment resolution in bp
#' @param parallel Whether to run samples in parallel
#' @param chr What chromosome to work on (leave as NULL to run on all chromosomes)
#' @param targets Samples/cells to shrink towards
#' @param cores How many cores to use when running samples in parallel
#' @param bootstrap Whether we should perform bootstrapping of inferred compartments
#' @param num.bootstraps How many bootstraps to run
#' @param genome What genome to work on ("hg19", "hg38", "mm9", "mm10")
#' @param other Another arbitrary genome to compute compartments on
#' @param group Whether to treat this as a group set of samples
#' @param boot.parallel Whether to run the bootstrapping in parallel
#' @param boot.cores How many cores to use for the bootstrapping
#'
#' @return A RaggedExperiment of inferred compartments
#' @import SummarizedExperiment
#' @import RaggedExperiment
#' @importFrom parallel mclapply
#' @importFrom methods as
#' @export
#'
#' @aliases getRNAABsignal
#'
#' @examples
#' if (requireNamespace("csaw", quietly = TRUE)) {
#'   data("k562_scatac_chr14", package = "compartmap")
#'   atac_compartments <- getATACABsignal(
#'     k562_scatac_chr14,
#'     parallel = FALSE,
#'     chr = "chr14",
#'     bootstrap = FALSE,
#'     genome = "hg19",
#'     group = TRUE
#'   )
#' }
getATACABsignal <- function(
  obj,
  res = 1e6,
  parallel = FALSE,
  chr = NULL,
  targets = NULL,
  cores = 2,
  bootstrap = TRUE,
  num.bootstraps = 100,
  genome = c("hg19", "hg38", "mm9", "mm10"),
  other = NULL,
  group = FALSE,
  boot.parallel = FALSE,
  boot.cores = 2
) {
  getCompartments(
    obj = obj,
    assay = "atac",
    res = res,
    parallel = parallel,
    chr = chr,
    targets = targets,
    cores = cores,
    bootstrap = bootstrap,
    num.bootstraps = num.bootstraps,
    boot.parallel = boot.parallel,
    boot.cores = boot.cores,
    genome = genome,
    group = group
  )
}

#' @describeIn getATACABsignal Alias for getATACABsignal
#'
getRNAABsignal <- getATACABsignal

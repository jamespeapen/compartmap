#' @title Estimate A/B compartments from single-cell RNA or ATAC sequencing data
#'
#' @description
#' \code{scCompartments} returns estimated A/B compartments from sc-seq data.
#'
#' @param obj Input SummarizedExperiment object
#' @param res Compartment resolution in bp
#' @param chr What chromosome to work on (leave as NULL to run on all chromosomes)
#' @param group Whether to treat this as a group set of samples
#' @param targets Samples/cells to shrink towards
#' @param bootstrap Whether we should perform bootstrapping of inferred compartments
#' @param num.bootstraps How many bootstraps to run
#' @param genome What genome to work on ("hg19", "hg38", "mm9", "mm10")
#' @param assay What type of single-cell assay is the input data ("atac" or "rna")
#' @param parallel Whether to run samples in parallel
#' @param cores How many cores to use when running samples in parallel
#' @param boot.parallel Whether to run the bootstrapping in parallel
#' @param boot.cores How many cores to use for the bootstrapping
#'
#' @return A RaggedExperiment of inferred compartments
#' @import SummarizedExperiment
#' @importFrom parallel mclapply
#' @import RaggedExperiment
#' @export
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' sc_compartments <- scCompartments(
#'   k562_scrna_chr14,
#'   chr = "chr14",
#'   parallel = FALSE,
#'   bootstrap = FALSE,
#'   genome = "hg19"
#' )
scCompartments <- function(
  obj,
  res = 1e6,
  chr = NULL,
  group = FALSE,
  targets = NULL,
  bootstrap = TRUE,
  num.bootstraps = 100,
  genome = c("hg19", "hg38", "mm9", "mm10"),
  assay = c("atac", "rna"),
  parallel = FALSE,
  cores = 2,
  boot.parallel = FALSE,
  boot.cores = 2
) {
  verifySE(obj)
  verifyCoords(obj)

  # which assay are we working on
  if (!all(assay %in% c("atac", "rna"))) stop("Supported assays are 'atac', and 'rna'.")
  assay <- tolower(match.arg(assay))
  verifyAssayNames(obj, assay = assay)
  getCompartments(
    obj = obj,
    assay = assay,
    res = res,
    chr = chr,
    targets = targets,
    parallel = parallel,
    cores = cores,
    bootstrap = bootstrap,
    num.bootstraps = num.bootstraps,
    boot.parallel = boot.parallel,
    boot.cores = boot.cores,
    genome = genome,
    group = group
  )
}

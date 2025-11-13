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
#' @param boot.parallel Whether to run the bootstrapping in parallel. See details.
#' @param BPPARAM BiocParallelParam object to use for parallelization. See details.
#'
#'
#' @return A RaggedExperiment of inferred compartments
#' @import SummarizedExperiment
#' @import RaggedExperiment
#' @importFrom BiocParallel bpparam
#' @export
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' sc_compartments <- scCompartments(
#'   k562_scrna_chr14,
#'   chr = "chr14",
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
  boot.parallel = FALSE,
  BPPARAM = bpparam()
) {
  verifySE(obj)
  verifyCoords(obj)

  bpparams <- get_nested_params(BPPARAM, boot.parallel)
  check_worker_count(bpparams)

  # which assay are we working on
  if (!all(assay %in% c("atac", "rna"))) stop("Supported assays are 'atac', and 'rna'.")
  assay <- tolower(match.arg(assay))
  verifyAssayNames(obj, assay = assay)
  getCompartments(
    obj = obj,
    res = res,
    chr = chr,
    group = group,
    targets = targets,
    bootstrap = bootstrap,
    num.bootstraps = num.bootstraps,
    genome = genome,
    assay = assay,
    boot.parallel = boot.parallel,
    bpparams = bpparams
  )
}

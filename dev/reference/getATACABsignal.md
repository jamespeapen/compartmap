# Estimate A/B compartments from ATAC-seq data

`getATACABsignal` returns estimated A/B compartments from ATAC-seq data.

## Usage

``` r
getATACABsignal(
  obj,
  res = 1000000,
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
)

getRNAABsignal(
  obj,
  res = 1000000,
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
)
```

## Arguments

- obj:

  Input SummarizedExperiment object

- res:

  Compartment resolution in bp

- parallel:

  Whether to run samples in parallel

- chr:

  What chromosome to work on (leave as NULL to run on all chromosomes)

- targets:

  Samples/cells to shrink towards

- cores:

  How many cores to use when running samples in parallel

- bootstrap:

  Whether we should perform bootstrapping of inferred compartments

- num.bootstraps:

  How many bootstraps to run

- genome:

  What genome to work on ("hg19", "hg38", "mm9", "mm10")

- other:

  Another arbitrary genome to compute compartments on

- group:

  Whether to treat this as a group set of samples

- boot.parallel:

  Whether to run the bootstrapping in parallel

- boot.cores:

  How many cores to use for the bootstrapping

## Value

A RaggedExperiment of inferred compartments

## Functions

- `getRNAABsignal()`: Alias for getATACABsignal

## Examples

``` r
if (requireNamespace("csaw", quietly = TRUE)) {
  data("k562_scatac_chr14", package = "compartmap")
  atac_compartments <- getATACABsignal(
    k562_scatac_chr14,
    parallel = FALSE,
    chr = "chr14",
    bootstrap = FALSE,
    genome = "hg19",
    group = TRUE
  )
}
#> Computing compartments for chr14
#> 108 bins created...
#> Calculating correlations...
#> Done...
#> Calculating eigenvectors.
#> Smoothing eigenvector.
#> Done smoothing.
```

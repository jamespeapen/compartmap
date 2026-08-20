#' Example scRNA-seq data for compartmap
#'
#' This object was generated using the K562 data from the STORM-seq paper
#' and pre-processed using the scran and scater packages and TF-IDF
#' transformed.
#'
#' @name k562_scrna_chr14
#' @docType data
#' @author Benjamin K Johnson \email{ben.johnson@vai.org}
#' @keywords data
#' @usage data(k562_scrna_chr14, package = "compartmap")
NULL

#' Example scRNA-seq data for compartmap
#'
#' This object was generated using the K562 data from the STORM-seq paper
#' and pre-processed using the scran and scater packages and are raw counts.
#'
#' @name k562_scrna_se_chr14
#' @docType data
#' @author Benjamin K Johnson \email{ben.johnson@vai.org}
#' @keywords data
#' @usage data(k562_scrna_raw, package = "compartmap")
NULL

#' Example scATAC-seq data for compartmap
#'
#' This data was generated using the data from the reference via bwa mem
#' and pre-processing the data using the csaw package.
#'
#' @name k562_scatac_chr14
#' @docType data
#' @author Benjamin K Johnson \email{ben.johnson@vai.org}
#' @references \url{https://www.ncbi.nlm.nih.gov//geo/query/acc.cgi?acc=GSE99172}
#' @keywords data
#' @usage data(k562_scatac_chr14, package = "compartmap")
NULL

#' Example SMART-seq3 scRNA-seq data for compartmap
#'
#' This object was generated using the HEK293T data from the SMART-seq3 paper
#'
#' @name ss3_umi_sce
#' @docType data
#' @author Benjamin K Johnson \email{ben.johnson@vai.org}
#' @references \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8735/}
#' @description Only keep chromosome 22 for the example
#' @keywords data
#' @usage data(ss3_umi_sce, package = "compartmap")
NULL

#' Example Illumina 450k methylation array data for compartmap
#'
#' This data was generated using the data from the reference via the
#' sesamize function from the SeSAMe package.
#'
#' @name array.data.chr14
#' @docType data
#' @author Benjamin K Johnson \email{ben.johnson@vai.org}
#' @references \url{https://f1000research.com/articles/5-1281/v3}
#' @keywords data
#' @usage data(array_data_chr14, package = "compartmap")
NULL

#' Human and mouse seqlengths as GRanges objects
#'
#' These objects were generated using the
#' BSgenome.Hsapiens.UCSC.hg19/BSgenome.Hsapiens.UCSC.hg38 and Mus.musculus
#' packages. The script used is found in the inst/scripts directory
#' @name seqlengths
#' @author Benjamin K Johnson \email{ben.johnson@vai.org}
#' @docType data
#' @keywords reference_data
#' @usage data(hg38.gr, package = "compartmap")
NULL

#' @name hg38.gr
#' @rdname seqlengths
"hg38.gr"

#' @name hg19.gr
#' @rdname seqlengths
"hg19.gr"

#' @name mm9.gr
#' @rdname seqlengths
"mm9.gr"

#' @name mm10.gr
#' @rdname seqlengths
"mm10.gr"

#' Human and mouse open sea CpGs as GRanges objects
#'
#' These objects were generated using the
#' BSgenome.Hsapiens.UCSC.hg19/BSgenome.Hsapiens.UCSC.hg38 and Mus.musculus
#' packages.
#'
#' @name openSeas
#' @author Benjamin K Johnson \email{ben.johnson@vai.org}
#' @docType data
#' @keywords reference_data
#' @usage data(openSeas.hg38, package = "compartmap")
NULL

#' @name openSeas.hg38
#' @rdname openSeas
"openSeas.hg38"

#' @name openSeas.hg19
#' @rdname openSeas
"openSeas.hg19"

#' @name openSeas.mm9
#' @rdname openSeas
"openSeas.mm9"

#' @name openSeas.mm10
#' @rdname openSeas
"openSeas.mm10"


#' Human and mouse genes as GRanges objects
#'
#' This object was generated using the
#' TxDb.Hsapiens.UCSC.hg19.knownGene/TxDb.Hsapiens.UCSC.hg38.knownGene/TxDb.Mmusculus.UCSC.mm9.knownGene/TxDb.Mmusculus.UCSC.mm10.knownGene
#' package. The script used for this object is found in the inst/scripts
#' directory
#'
#' @name genes
#' @author James Eapen \email{james.eapen@vai.org}
#' @docType data
#' @keywords reference_data
#' @usage data(hg38.tx.gr, package = "compartmap")
NULL

#' @name hg38.tx.gr
#' @rdname genes
"hg38.tx.gr"

#' @name hg19.tx.gr
#' @rdname genes
"hg19.tx.gr"

#' @name mm9.tx.gr
#' @rdname genes
"mm9.tx.gr"

#' @name mm10.tx.gr
#' @rdname genes
"mm10.tx.gr"

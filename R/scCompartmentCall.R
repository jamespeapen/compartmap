#' @rdname CompartmentCall
#'
#' @param re A `RaggedExperiment` of single-cell compartment calls
#' @param unitarize Whether to unitarize the singular values for each of the inputs calls
#'
#' @importFrom data.table melt as.data.table
#' @keywords CompartmentCall
#' @export
scCompartmapCall <- new_class(
  "scCompartmapCall",
  parent = MultiCompartmapCall,
  properties = list(
    colnames = class_character,
    mat = new_S3_class(c("matrix", "array"))
  ),
  constructor = function(re, res, name, unitarized = FALSE, unitarize = FALSE) {
    grlist <- condenseSE(re)
    gen <- genome(re)
    pcs <- lapply(grlist, function(i) {
      mcols(i)[, 'pc']
    })
    mat <- do.call(cbind, pcs)

    if (unitarize && !unitarized) {
      mat <- apply(mat, 2, .unitarize)
    } else if (unitarize && unitarized) {
      message("Already unitarized")
    }

    df <- melt(
      as.data.table(mat)[, n := .I],
      id.vars = "n",
      variable.name = "name",
      value.name = "pc"
    )

    gr <- GRanges(rownames(mat))
    genome(gr) <- gen

    new_object(
      S7_object(),
      name = name,
      gr = gr,
      df = df,
      res = res,
      unitarized = unitarized,
      seqinfo = methods::selectMethod('seqinfo', "GRanges")(gr),
      colnames = colnames(mat),
      mat = mat,
      filtered = FALSE,
      filter_threshold = 0
    )
  }
)
S4_register(scCompartmapCall)

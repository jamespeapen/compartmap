#' Compartment Call objects for analysis (experimental)
#'
#' These S7 classes hold compartment call values, and their genomic region and
#' resolution metadata as well as methods to aid analysis and visualization.
#' `CompartmentCall()` is primarily a parent class from which
#' `CompartmapCall()` and `MultiCompartmapCall()` are derived.
#'
#' @details
#'
#' ## Objects
#'
#' `CompartmentCall` is constucted from a vector of the call values, a
#' `GenomicRanges::GRanges()` object of the bins, and the resolution of the
#' bins. Given a `GRanges` of bin coordinates, it can be used to store singular
#' values/eigenvectors from Hi-C compartment analysis. Along with
#' `CompartmapCall` objects, they can be combined in a `MultiCompartmapCall` to
#' compare against each other. All `CompartmentCall` methods work on every
#' other *compartmap* S7 object.
#'
#' `CompartmapCall()` is constructed from the `GRanges` output of
#' `scCompartments(group = TRUE)`. `MultiCompartmapCall()` holds multiple
#' `CompartmapCall()` objects at the same resolution, and is constructed with a
#' list of `CompartmentCall` or `CompartmapCall` objects.
#'
#' For single-cell level compartment inferences use `scCompartmapCall()`, which
#' inherits from `MultiCompartmapCall`. This class holds multiple single-cell
#' level compartment calls, taking the `RaggedExperiment::RaggedExperiment()`
#' from `scCompartments()` as its constructor. Since it inherits from
#' `MultiCompartmapCall`, all `MultiCompartmapCall` methods work on
#' `scCompartmapCall` objects.
#'
#' All objects may be subset to required indices with `object[row indices,
#' column indices]`. Their dimensions can be accessed with `dim()`, `nrow()`
#' and `ncol()`.
#'
#' ## Properties and accessors
#'
#' These properties may be accessed with `@` like S4 slots, or using their
#' accessor functions. Some functions share methods with other Bioconductor
#' classes like `GRanges` and `SummarizedExperiment`.
#'
#' - `name`, `get_name()`: The name of the object, used to identify it, useful
#' when making `MultiCompartmapCall` objects
#' - `gr`, `granges()`: A `GRanges` object of the compartment bins
#' - `df`, `DF()`: a `data.table` of bin indices in column `n` and compartment
#' call singular values in column `cscore`. For `MultiCompartmapCall` and
#' `scCompartmapCall` objects, this is in a tidy format, with an additional
#' `name` column.
#' - `res`, `resolution()`: The genomic resolution at which the compartments
#' were called
#' - `unitarized`, `is_unitarized()`: Whether the singular values have been
#' unitarized
#' - `seqinfo`, `seqinfo()`: Seqinfo for the object's `GRanges`
#' - `mat`, `mat()`: A matrix of genomic bins by singular values - a wide
#' matrix format of the `df` slot. This is used to calculate correlation and
#' agreement between cells in a `scCompartmapCall` and groups in a
#' `MultiCompartmapCall` object.
#'
#' @param cscore The singular values from a compartment call
#' @param res The binning resolution used
#' @param gr The GRanges of the bins or the output of `scCompartments` or
#' `getArrayCompartments` containing the 'cscore' column
#' @param assay What assay is this from: RNA, ATAC, methylation, Hi-C?
#' @param name An identifier for the object. For `CompartmentCall` and
#' `CompartmapCall`, this is becomes the column name for the object when added
#' to a `MultiCompartmapCall` object
#' @param unitarized Whether the singular values have been unitarized
#'
#' @import S7
#' @importFrom data.table data.table .I .N :=
#' @keywords CompartmentCall
#' @export
CompartmentCall <- new_class(
  "CompartmentCall",
  properties = list(
    name = class_character,
    gr = methods::getClass("GRanges"),
    df = methods::getClass("data.table"),
    res = class_numeric,
    assay = class_character,
    unitarized = class_logical,
    filtered = class_logical,
    filter_threshold = class_numeric,
    seqinfo = methods::getClass("Seqinfo")
  ),
  constructor = function(cscore, res, gr, assay, name = NULL, unitarized = FALSE) {
    df <- data.table(cscore = as.vector(cscore))[, .(n = .I, pos = start(gr), cscore, name = name)]
    new_object(
      S7_object(),
      name = name %||% shQuote(substitute(gr), "cmd2"),
      gr = granges(gr),
      df = df,
      res = res,
      assay = assay,
      unitarized = unitarized,
      filtered = FALSE,
      filter_threshold = 0,
      seqinfo = methods::selectMethod('seqinfo', "GRanges")(gr)
    )
  },
  validator = function(self) {
    if (anyNA(genome(self@gr))) {
      "`gr`'s `genome` must be specified. Set with `genome(gr)` <- [genome]"
    }
  }
)
S4_register(CompartmentCall)

# GETTERS ==================================================================={{{

#' Get the `@df` slot from a CompartmentCall object.
#'
#' `n`: bin indices corresponding to indices of the `GRanges` object in `@gr`
#' `cscore`: compartment call singular values
#' `name`: The name of the individual `CompartmentCall` in a
#' `MultiCompartmapCall` object
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
DF <- new_generic("DF", "x", function(x) S7_dispatch())
method(DF, CompartmentCall) <- function(x) x@df[]

#' Get dimensions of `CompartmapCall` matrix
#' @concept s7getters
#' @rdname dim
#' @export
method(dim, CompartmentCall) <- function(x) c(length(x@gr), 1)


#' Get GRanges of the `CompartmentCall`
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
method(granges, CompartmentCall) <- function(x) x@gr

#' Get Seqinfo of the `CompartmentCall`
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
method(seqinfo, CompartmentCall) <- function(x) x@seqinfo

#' Get 'seqlevels' of the `CompartmentCall`
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
method(seqlevels, CompartmentCall) <- function(x) {
  methods::selectMethod('seqlevels', 'GRanges')(x)
}

#' Get the resolution of the `CompartmentCall`
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
get_resolution <- new_generic("get_resolution", "x", function(x) S7_dispatch())
method(get_resolution, CompartmentCall) <- function(x) x@res

#' Get the assay of the `CompartmentCall`
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
get_assay <- new_generic("get_assay", "x", function(x) S7_dispatch())
method(get_assay, CompartmentCall) <- function(x) x@assay

#' Check if the `CompartmentCall` was unitarized
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
is_unitarized <- new_generic("is_unitarized", "x", function(x) S7_dispatch())
method(is_unitarized, CompartmentCall) <- function(x) x@unitarized

#' Get the name of the `CompartmentCall` object
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
get_name <- new_generic("get_name", "x", function(x) S7_dispatch())
method(get_name, CompartmentCall) <- function(x) x@name

#' Check if the `CompartmentCall` was filtered.
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
is_filtered <- new_generic("is_filtered", "x", function(x) S7_dispatch())
method(is_filtered, CompartmentCall) <- function(x) x@filtered

#' Check the threshold at which an object's calls were filtered, returning 0 if not filtered
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
get_filter_threshold <- new_generic("filter_threshold", "x", function(x) S7_dispatch())
method(get_filter_threshold, CompartmentCall) <- function(x) x@filter_threshold

#' Find indices of open and closed compartments
#'
#' @param x A `CompartmentCall` object
#'
#' @return
#' - `is_open`: boolean vector where TRUE means open
#' - `is_closed`: boolean vector where TRUE means closed
#' @concept s7getters
#' @export
is_open <- new_generic("is_open", "x", function(x) S7_dispatch())
method(is_open, CompartmentCall) <- function(x) x@df[, cscore > 0]

#' @rdname is_open
#' @export
is_closed <- new_generic("is_closed", "x", function(x) S7_dispatch())
method(is_closed, CompartmentCall) <- function(x) Negate(is_open)(x)

# }}}

# SUBSETTERS ==============================================================={{{

#' Subset rows of a `CompartmentCall` object
#'
#' @param x A `CompartmentCall` object
#' @param i Row indices to subset
#'
#' @concept s7ranges
#' @export
#' @keywords internal
`[.compartmap::CompartmentCall` <- function(x, i = NULL) {
  i <- i %||% seq_len(nrow(x@mat))
  x@df <- x@df[i]
  x@df[, n := .I][]
  x@gr <- x@gr[i]
  Seqinfo::seqlevels(x@gr) <- Seqinfo::seqlevelsInUse(x@gr)
  x
}

#' Subset the `CompartmentCall` object by chromosome
#'
#' @param x A `CompartmentCall` object
#' @param chr A string vector of chromosomes to subset to
#'
#' @concept s7ranges
#' @export
subset_chr <- new_generic("subset_chr", "x", function(x, chr) S7_dispatch())
method(subset_chr, CompartmentCall) <- function(x, chr) {
  ind <- as.vector(seqnames(x@gr) %gin% chr)
  x[which(ind)]
}

#' Filter to bins with call values greater than or equal to a threshold value
#'
#' @param x A `CompartmentCall` object
#' @param threshold The absolute value to use for filtering. Rows where any
#' value is less than this threshold are dropped
#'
#' @concept s7analysis
#' @export
filter <- new_generic("filter", "x", function(x, threshold = 0.02) S7_dispatch())
method(filter, CompartmentCall) <- function(x, threshold = 0.02) {
  filter_rows <- x@df[, abs(cscore) >= threshold]
  x <- x[filter_rows]
  x@filtered <- TRUE
  x@filter_threshold <- threshold
  x
}

# }}}

#' Find overlaps between `CompartmentCall` objects
#'
#' @param subject A `CompartmentCall` object
#' @param query A string vector of chromosomes to subset to
#'
#' @concept s7ranges
#' @export
method(findOverlaps, list(CompartmentCall, CompartmentCall)) <- function(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
) {
  query <- query@gr
  subject <- subject@gr
  methods::callGeneric()
}

#' Unitarize the singular values in a `CompartmentCall` or `MultiCompartmapCall`
#'
#' @param x A `CompartmentCall` object
#' @param medianCenter Whether to center the singular values on their median
#'
#' @concept s7analysis
#' @export
unitarize <- new_generic("unitarize", "x", function(x, medianCenter = TRUE) S7_dispatch())
method(unitarize, CompartmentCall) <- function(x, medianCenter = TRUE) {
  if (is_unitarized(x)) {
    message(get_name(x), " is already unitarized")
    return(x)
  }

  df <- x@df
  if (inherits(x, "compartmap::MultiCompartmapCall")) {
    x@mat <- apply(x@mat, 2, .unitarize, medianCenter = medianCenter)
    x@df <- df[,
      cscore := .unitarize(cscore, medianCenter = medianCenter),
      by = name
    ][]
  } else {
    x@df <- df[, cscore := .unitarize(cscore, medianCenter = medianCenter)]
  }

  x@unitarized <- TRUE
  x
}

# SIGN CORRECTION ==========================================================={{{

#' Flip the singular values signs in a `CompartmentCall`
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7analysis
#' @export
flip <- new_generic("flip", "x", function(x) S7_dispatch())
method(flip, CompartmentCall) <- function(x) {
  x@df <- data.table::copy(x@df)[, cscore := -cscore]
  x
}

#' Correct the sign of the compartment call vector
#'
#' @param x A `CompartmentCall` object
#' @param na.rm Whether to remove NAs. The presence of NAs can reduce the
#' accuracy of the sign flip as the open vs closed gene densities cannot be
#' computed correctly
#'
#' @concept s7analysis
#' @export
fix_sign <- new_generic("fix_sign", "x", function(x, na.rm = FALSE) S7_dispatch())
method(fix_sign, CompartmentCall) <- function(x, na.rm = FALSE) {
  nas_present <- anyNA(x@df)
  if (!na.rm && nas_present) {
    warning(
      "NAs found - flipping may fail. If you think the NAs will not affect the distribution of expression values set `na.rm = TRUE`"
    )
  }
  gr <- x@gr
  gr$cscore <- x@df[, cscore]
  gr_full <- gr
  if (na.rm) {
    gr_full <- gr[!is.na(gr$cscore)]
  }
  if (flipSign(gr_full, genome(gr), x@assay)) {
    x <- flip(x)
  }
  x
}

# }}}

#' Fill missing genomic bins in `CompartmentCalls` using a reference GRanges
#'
#' Compartmap may drop genomic bins with insufficient data and the resulting
#' GRanges object may not have all the bins of the region it was run which
#' means using the same region and resolution on different inputs does not
#' guarantee the same output bins. Having different set of bins between calls
#' despite the same input regions and resolution prevents the creation of
#' `MultiCompartmapCall` objcets that expect the same GRanges across all input
#' `CompartmentCall` objects. `fill_missing()` adds the missing bins according
#' to a larger reference set of bins, filling missing data with NA. All
#' bins in `x` must be present in the reference bins.
#'
#' @param x A `CompartmentCall` object
#' @param ref.gr The `GRanges` object to use as the full reference set of
#' regions
#'
#' @concept s7ranges
#' @export
fill_missing <- new_generic("fill_missing", "x", function(x, ref.gr) S7_dispatch())
method(fill_missing, CompartmentCall) <- function(x, ref.gr) {
  ref_length <- length(ref.gr)
  stopifnot("Reference GRanges is not bigger than CompartmentCall object" = ref_length >= length(x@gr))
  stopifnot("All CompartmentCall bins must be present in the reference GRanges" = all(x@gr %gin% ref.gr))

  ol <- findOverlaps(x@gr, ref.gr)
  mcols(ref.gr)[, c("cscore", "name")] <- NA
  ref.gr[subjectHits(ol)]$cscore <- x@df[queryHits(ol), cscore]
  ref.gr[subjectHits(ol)]$name <- x@df[queryHits(ol), name]
  x@gr <- granges(ref.gr)
  x@df <- as.data.table(mcols(ref.gr))[, .(n = .I, pos = start(ref.gr), cscore, name)]
  x
}


#' Get the difference between two CompartmentCall objects call values
#' @param x,y CompartmentCall objects to compare
#' @concept s7analysis
#' @rdname CompartmentCall_difference
#' @export
`-.compartmap::CompartmentCall` <- function(x, y) {
  stopifnot("Both objects must have the same GRanges" = length(x@gr) == length(y@gr))
  stopifnot("Both objects must have the same @df" = nrow(x@df) == nrow(y@df))
  df <- x@df
  cscore1 <- DF(x)[, cscore]
  cscore2 <- DF(y)[, cscore]
  x@df <- x@df[, `:=`(cscore = cscore1 - cscore2, name = paste(x@name, "-", y@name))]
  x
}

#' Plot singular values of a `CompartmentCall` object
#'
#' @param x `CompartmentCall` or `CompartmapCall` object
#' @param ... Placeholder for the `plot` generic - arguments have not effect
#' @param type Whether to plot the singular values as `"line"` or `"bar"` plots.
#' @param res The resolution to round the genomic coordinates to (kilobase:
#' "kb" or megabase: "mb")
#' @param width The width of the `geom_line` if `type = "line"` or the width
#' of the bar if `type = "bar"` in the plot
#' @param ylim Upper and lower bound for the y-axis
#'
#' @importFrom ggplot2 ggplot geom_line geom_bar scale_x_continuous scale_y_continuous theme element_text
#' @concept plotting
#' @export
`plot.compartmap::CompartmentCall` <- function(
  x,
  ...,
  type = "line",
  res = "kb",
  width = 0.5,
  ylim = NULL
) {
  . <- NULL
  cscore <- NULL

  pd <- x@df
  lim <- ylim %||% range(pd$cscore) |> abs() |> max() * c(-1, 1)

  p <- switch(
    type,
    line = {
      ggplot(pd, aes(x = pos, y = cscore, color = name, group = name)) +
        geom_line()
    },
    bar = {
      ggplot(pd, aes(x = pos, y = cscore, group = name, fill = cscore > 0)) +
        geom_col() +
        facet_grid(rows = vars(name)) +
        scale_fill_manual(values = c("deeppink4", "grey50")) +
        theme(legend.position = "none")
    }
  )

  p +
    scale_y_continuous(limits = ylim) +
    scale_x_continuous(labels = \(x) x / 1e6) +
    labs(x = paste(gsub("chr", "Chromosome ", seqlevels(x)), "(Mb)"))
}

.resolution <- function(res) {
  scales::number(res, suffix = "b", scale_cut = scales::cut_short_scale())
}

# Store the BiocGenerics %in% to differentiate from base::`%in%`
`%gin%` <- function(a, b) {
  BiocGenerics::match(a, b, nomatch = 0L) > 0L
}

#' Print `CompartmentCall`
#' @export
#' @keywords internal
method(print, CompartmentCall) <- function(x, ...) {
  message(.print_CompartmentCall(x))
}

.print_CompartmentCall <- function(x) {
  properties <- props(x)
  class_type <- gsub("compartmap::", "", class(x)[1])

  sprintf(
    "<%s object>
  @name        : %s
  @assay       : %s
  @res         : %s
  @gr          : GRanges with %d bins
  @df          : data.table of compartment calls (n = bin index, cscore = singular values)
  @unitarized  : %s
  @filtered    : %s%s",
    class_type,
    x@name,
    x@assay,
    .resolution(x@res),
    length(x@gr),
    x@unitarized,
    x@filtered,
    ifelse(x@filtered, paste(", for absolute values >=", x@filter_threshold), "")
  )
}

#' @rdname CompartmentCall
#'
#' @keywords CompartmentCall
#' @export
CompartmapCall <- new_class(
  "CompartmapCall",
  parent = CompartmentCall,
  constructor = function(gr, res, assay, name = NULL, unitarized = FALSE) {
    if ("score" %in% colnames(mcols(gr))) {
      df <- data.table(cscore = gr$score)
    } else {
      df <- data.table(cscore = gr$cscore)
    }
    df <- df[, `:=`(n = .I, name = name, pos = start(gr))][, .(n, pos, cscore, name)]
    new_object(
      S7_object(),
      name = name %||% shQuote(substitute(gr), "cmd2"),
      gr = granges(gr),
      df = df,
      res = res,
      assay = assay,
      unitarized = unitarized,
      filtered = FALSE,
      filter_threshold = 0,
      seqinfo = methods::selectMethod("seqinfo", "GRanges")(gr)
    )
  }
)
S4_register(CompartmapCall)

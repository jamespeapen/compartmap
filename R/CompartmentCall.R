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
#' call singular values in column `pc`. For `MultiCompartmapCall` and
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
#' @param pc The singular values from a compartment call
#' @param res The binning resolution used
#' @param gr The GRanges of the bins or the output of `scCompartments` or
#' `getArrayCompartments` containing the 'pc' column
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
    unitarized = class_logical,
    filtered = class_logical,
    filter_threshold = class_numeric,
    seqinfo = methods::getClass("Seqinfo")
  ),
  constructor = function(pc, res, gr, name = NULL, unitarized = FALSE) {
    df <- data.table(pc = as.vector(pc))[, .(n = .I, pc, name = name)]
    new_object(
      S7_object(),
      name = name %||% shQuote(substitute(gr), "cmd2"),
      gr = granges(gr),
      df = df,
      res = res,
      unitarized = unitarized,
      filtered = FALSE,
      filter_threshold = 0,
      seqinfo = methods::selectMethod('seqinfo', "GRanges")(gr)
    )
  },
  validator = function(self) {
    if (is.na(genome(self@gr))) {
      "`gr`'s `genome` must be specified. Set with `genome(gr)` <- [genome]"
    }
  }
)
S4_register(CompartmentCall)

# GETTERS ==================================================================={{{

#' Get the `@df` slot from a CompartmentCall object.
#'
#' `n`: bin indices corresponding to indices of the `GRanges` object in `@gr`
#' `pc`: compartment call singular values
#' `name`: The name of the individual `CompartmentCall` in a
#' `MultiCompartmapCall` object
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
DF <- new_generic("DF", "x", function(x) {
  S7_dispatch()
})
method(DF, CompartmentCall) <- function(x) {
  x@df[]
}

#' @concept s7getters
`nrow.compartmap::CompartmentCall` <- function(x) {
  length(x@gr)
}

#' @concept s7getters
`ncol.compartmap::CompartmentCall` <- function(x) {
  1
}

#' @concept s7getters
#' Get GRanges of the `CompartmentCall`
#'
#' @param x A `CompartmentCall` object
#'
#' @export
method(granges, CompartmentCall) <- function(x) {
  x@gr
}

#' Get Seqinfo of the `CompartmentCall`
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
method(seqinfo, CompartmentCall) <- function(x) {
  x@seqinfo
}


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
resolution <- new_generic("resolution", "x", function(x) {
  S7_dispatch()
})
method(resolution, CompartmentCall) <- function(x) {
  x@res
}

#' Check if the `CompartmentCall` was unitarized
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
is_unitarized <- new_generic("is_unitarized", "x", function(x) {
  S7_dispatch()
})
method(is_unitarized, CompartmentCall) <- function(x) {
  x@unitarized
}

#' Get the name of the `CompartmentCall` object
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
get_name <- new_generic("get_name", "x", function(x) {
  S7_dispatch()
})
method(get_name, CompartmentCall) <- function(x) {
  x@name
}

#' Check if the `CompartmentCall` was filtered.
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
is_filtered <- new_generic("is_filtered", "x", function(x) {
  S7_dispatch()
})
method(is_filtered, CompartmentCall) <- function(x) {
  x@filtered
}

#' Check the threshold at which an object's calls were filtered, returning 0 if not filtered
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7getters
#' @export
get_filter_threshold <- new_generic("filter_threshold", "x", function(x) {
  S7_dispatch()
})
method(get_filter_threshold, CompartmentCall) <- function(x) {
  x@filter_threshold
}

# }}}

# SUBSETTERS ==============================================================={{{

#' Subset rows of a `CompartmentCall` object
#'
#' @param x A `CompartmentCall` object
#' @param i Row indices to subset
#'
#' @concept s7ranges
#' @export
`[.compartmap::CompartmentCall` <- function(x, i = NULL) {
  i <- i %||% seq_len(nrow(x@mat))
  x@df <- x@df[i]
  x@df[, n := .I][]
  x@gr <- x@gr[i]
  seqlevels(x@gr) <- seqlevelsInUse(x@gr)
  x
}

#' Subset the `CompartmentCall` object by chromosome
#'
#' @param x A `CompartmentCall` object
#' @param chr A string vector of chromosomes to subset to
#'
#' @concept s7ranges
#' @export
subset_chr <- new_generic("subset_chr", "x", function(x, chr) {
  S7_dispatch()
})
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
filter <- new_generic("filter", "x", function(x, threshold = 0.02) {
  S7_dispatch()
})
method(filter, CompartmentCall) <- function(x, threshold = 0.02) {
  filter_rows <- x@df[, abs(pc) >= threshold]
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
unitarize <- new_generic("unitarize", "x", function(x, medianCenter = TRUE) {
  S7_dispatch()
})
method(unitarize, CompartmentCall) <- function(x, medianCenter = TRUE) {
  if (is_unitarized(x)) {
    message(get_name(x), " is already unitarized")
    return(x)
  }

  df <- x@df
  if (inherits(x, "compartmap::MultiCompartmapCall")) {
    x@mat <- apply(x@mat, 2, .unitarize, medianCenter = medianCenter)
    x@df <- x@df[, .(n, pc = .unitarize(pc, medianCenter = medianCenter)), by = name][, .(n, pc, name)]
  } else {
    x@df <- x@df[, .(n, pc = .unitarize(pc, medianCenter = medianCenter), name)]
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
flip <- new_generic("flip", "x", function(x) {
  S7_dispatch()
})
method(flip, CompartmentCall) <- function(x) {
  x@df <- x@df[, .(n, pc = -pc, name)]
  x
}

#' Correct the sign of the compartment call vector
#'
#' @param x A `CompartmentCall` object
#'
#' @concept s7analysis
#' @export
fix_sign <- new_generic("fix_sign", "x", function(x) {
  S7_dispatch()
})
method(fix_sign, CompartmentCall) <- function(x) {
  gr <- x@gr
  gr$pc <- x@df[, pc]
  if (flipSign(gr, genome(gr))) {
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
fill_missing <- new_generic("flip", "x", function(x, ref.gr) {
  S7_dispatch()
})
method(fill_missing, CompartmentCall) <- function(x, ref.gr) {
  ref_length <- length(ref.gr)
  stopifnot("Reference GRanges is not bigger than CompartmentCall object" = ref_length > length(x@gr))
  stopifnot("All CompartmentCall bins must be present in the reference GRanges" = all(x@gr %gin% ref.gr))

  ref_idx <- seq_len(ref_length)
  df <- data.table(n = ifelse(ref.gr %gin% x@gr, ref_idx, NA))
  df[!is.na(n), `:=`(pc = x@df$pc, name = x@name)]
  df[, n := .I][]
  x@gr <- ref.gr
  x@df <- df
  x
}


#' Get the difference between two CompartmentCall objects call values
#' @param x, y CompartmentCall objects to compare
#' @concept s7analysis
#' @rdname CompartmentCall_difference
#' @export
`-.compartmap::CompartmentCall` <- function(x, y) {
  stopifnot("Both objects must have the same GRanges" = length(x@gr) == length(y@gr))
  stopifnot("Both objects must have the same @df" = nrow(x@df) == nrow(y@df))
  df <- x@df
  pc1 <- DF(x)[, pc]
  pc2 <- DF(y)[, pc]
  x@df <- x@df[, .(n, pc = pc1 - pc2, name = paste(x@name, "-", y@name))]
  x
}

#' Plot singular values of a `CompartmentCall` object
#'
#' @param x `CompartmentCall` or `CompartmapCall` object
#' @param ... Placeholder for the `plot` generic - arguments have not effect
#' @param type Whether to plot the singular values as `"line"` or `"bar"` plots.
#' @param label_coords Label the x-axis with genomic coordinates. Uses a
#' numeric index when set to `FALSE`. Using coordinate labels can severely
#' crowd the x-axis, especially with Kb-resolution calls.
#' @param res The resolution to round the genomic coordinates to (kilobase:
#' "kb" or megabase: "mb")
#' @param width The width of the `geom_line` if `type = "line"` or the width
#' of the bar if `type = "bar"` in the plot
#' @param ylim Upper and lower bound for the y-axis
#'
#' @importFrom ggplot2 ggplot geom_line geom_bar scale_y_continuous theme element_text
#' @concept plotting
#' @export
`plot.compartmap::CompartmentCall` <- function(
  x,
  ...,
  type = "line",
  label_coords = FALSE,
  res = "kb",
  width = 0.5,
  ylim = c(-1, 1)
) {
  . <- NULL
  pc <- NULL
  if (label_coords) {
    coord_pd <- x@df[, .(n, pc, coord = grscale(x@gr, res))]
    p <- ggplot(coord_pd, aes(x = coord, y = pc, group = 1)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  } else {
    p <- ggplot(x@df, aes(x = n, y = pc))
  }
  p <- p + scale_y_continuous(limits = ylim)

  switch(
    type,
    line = p + geom_line(linewidth = width),
    bar = p + geom_bar(stat = "identity", width = width)
  )
}

grscale <- function(gr, res) {
  scale_factor <- switch(tolower(res), kb = list(1e5, "Kb"), mb = list(1e6, "Mb"))
  start_scaled <- start(gr) / as.numeric(scale_factor[1])
  end_scaled <- round(end(gr) / as.numeric(scale_factor[1], 4))
  paste0(seqlevels(gr), ":", start_scaled, "-", end_scaled, " ", scale_factor[2])
}

.resolution <- function(res) {
  fct <- res / 1e5
  if ((any(fct < 1) | any(fct >= 100))) {
    stop("Unsupported resolution")
  }
  ifelse(fct < 10, paste(fct, "Kb"), paste(fct / 10, "Mb"))
}

# Store the BiocGenerics %in% to differentiate from base::`%in%`
`%gin%` <- function(a, b) {
  BiocGenerics::match(a, b, nomatch = 0L) > 0L
}

#' Print `CompartmentCall`
method(print, CompartmentCall) <- function(x, ...) {
  message(.print_CompartmentCall(x))
}

.print_CompartmentCall <- function(x) {
  properties <- props(x)
  class_type <- gsub("compartmap::", "", class(x)[1])

  sprintf(
    "<%s object>
  @name        : %s
  @res         : %s
  @gr          : GRanges with %d bins
  @df          : data.table of compartment calls (n = bin index, pc = singular values)
  @unitarized  : %s
  @filtered    : %s%s",
    class_type,
    x@name,
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
  constructor = function(gr, res, name = NULL, unitarized = FALSE) {
    df <- data.table(pc = gr$pc)[, `:=`(n = .I, name = name)][, .(n, pc, name)]
    new_object(
      S7_object(),
      name = name %||% shQuote(substitute(gr), "cmd2"),
      gr = granges(gr),
      df = df,
      res = res,
      unitarized = unitarized,
      filtered = FALSE,
      filter_threshold = 0,
      seqinfo = methods::selectMethod("seqinfo", "GRanges")(gr)
    )
  }
)
S4_register(CompartmapCall)

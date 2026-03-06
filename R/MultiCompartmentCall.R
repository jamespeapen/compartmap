#' @rdname CompartmentCall
#'
#' @param ccalls A list of `CompartmapCalls` to combine
#' @param name An identifier for this set of compartment calls
#' @param unitarize Whether to unitarize the singular values for each of the inputs calls
#' @param filter Whether to filter the singular values
#' @param filter_threshold The threshold at which to filter the singular values
#'
#' @importFrom data.table rbindlist dcast
#' @keywords CompartmentCall
#' @export
MultiCompartmapCall <- new_class(
  "MultiCompartmapCall",
  parent = CompartmentCall,
  properties = list(
    colnames = class_character,
    mat = new_S3_class(c("matrix", "array"))
  ),
  constructor = function(ccalls, name, unitarize = FALSE, filter = FALSE, filter_threshold = 0.02) {
    all_same <- function(prop_get, err) {
      unique_property <- unique(unlist(lapply(ccalls, prop_get)))
      if (length(unique_property) != 1) {
        stop(err)
      }
      unique_property
    }
    unitarize_all <- function(i) {
      if (!is_unitarized(i)) {
        return(unitarize(i))
      }
      i
    }

    unique_res <- all_same(resolution, "All resolutions must be the same")
    unique_gr <- all_same(granges, "All GRanges must contain the same ranges")
    all_filtered <- all_same(is_filtered, "All objects must be filtered")
    all_ft <- all_same(get_filter_threshold, "All objects must be filtered to the same threshold")

    all_unitarized <- tryCatch(
      all_same(
        is_unitarized,
        "All calls must be either unitarized or non-unitarized - set 'unitarize = TRUE' to unitarize them"
      ),
      error = function(e) {
        if (unitarize) {
          FALSE
        } else {
          stop(e)
        }
      }
    )

    if (unitarize) {
      ccalls <- lapply(ccalls, unitarize)
      all_unitarized <- TRUE
      all_ft <- filter_threshold
    }

    colorder <- sapply(ccalls, get_name)
    df <- rbindlist(lapply(ccalls, function(i) {
      DF(i)[, name := get_name(i)][, name := factor(name, levels = colorder)][]
    }))
    mat <- as.matrix(dcast(df, n ~ name, value.var = "pc")[, -1])

    obj <- new_object(
      S7_object(),
      name = name,
      gr = unique_gr[[1]],
      df = df,
      res = unique_res,
      unitarized = all_unitarized,
      filtered = all_filtered,
      filter_threshold = all_ft,
      seqinfo = methods::selectMethod('seqinfo', "GRanges")(unique_gr[[1]]),
      colnames = unlist(lapply(ccalls, get_name)),
      mat = mat
    )

    if (filter) {
      obj <- filter(obj, threshold = filter_threshold)
    }
    obj
  },
  validator = function(self) {
    if (length(self@colnames) != length(unique(self@colnames))) {
      "All names must be unique"
    }
  }
)
S4_register(MultiCompartmapCall)

# GETTERS ==================================================================={{{
#' Get dimensions of `MultiCompartmapCall` matrix
#' @concept s7getters
#' @export
method(dim, MultiCompartmapCall) <- function(x) dim(x@mat)


#' Get column names of `MultiCompartmapCall` matrix
#' @concept s7getters
#' @export
method(names, MultiCompartmapCall) <- function(x) colnames(x@mat)


#' Get the data matrix from a `MultiCompartmapCall` or `scCompartmapCall`
#'
#' @param x A `MultiCompartmapCall` object
#'
#' @concept s7getters
#' @export
mat <- new_generic("mat", "x", function(x) S7_dispatch())
method(mat, MultiCompartmapCall) <- function(x) x@mat
# }}}

# SUBSETTERS ==============================================================={{{

#' Subset `MultiCompartmapCall` object
#'
#' @param x A `CompartmentCall` object
#' @param i Rows indices to subset
#' @param j Column indices or names to subset
#'
#' @concept s7ranges
#' @export
`[.compartmap::MultiCompartmapCall` <- function(x, i = NULL, j = NULL) {
  i <- i %||% seq_len(nrow(x@mat))
  j <- j %||% seq_len(ncol(x@mat))

  if (is.logical(i)) {
    i <- which(i)
  } else if (!is.numeric(i)) {
    stop("Unsupported type for subsetting columns")
  }
  stopifnot("i exceeds the row count of x" = i <= nrow(x@mat))

  if (is.numeric(j)) {
    if (length(j) > ncol(x@mat)) {
      stop("j exceeds the column count of x")
    }
    subset_names <- x@df[, unique(name)][j]
  } else if (is.character(j)) {
    stopifnot("One or more columns not found in the data" = all(j %in% colnames(x@mat)))
    subset_names <- j
  } else {
    stop("Unsupported type for subsetting columns")
  }

  x@gr <- x@gr[i]
  seqlevels(x@gr) <- seqlevelsInUse(x@gr)
  x@mat <- x@mat[i, j, drop = FALSE]
  x@df <- x@df[n %in% i & name %in% subset_names]
  x@df[, n := seq_len(.N), by = name][]
  x
}

# }}}

#' Print a `CompartmentCall` object
#'
#' @param x A `CompartmentCall` object
#'
#' @export
method(print, MultiCompartmapCall) <- function(x, ...) {
  column_label <- ifelse(inherits(x, "compartmap::scCompartmapCall"), "cells", "samples")
  msg <- message(
    .print_CompartmentCall(x),
    sprintf("\n  @mat         : %d bins x %d %s", nrow(x@mat), ncol(x@mat), column_label)
  )
}

#' @concept s7analysis
#' @export
method(filter, MultiCompartmapCall) <- function(x, threshold = 0.02) {
  filter_rows <- apply(x@mat, 1, function(i) {
    all(abs(i) >= threshold)
  })
  x <- x[filter_rows]
  x@filtered <- TRUE
  x@filter_threshold <- threshold
  x
}

#' Compute agreement between compartment calls
#'
#' Calculates the proportion of calls with the same sign for every pair of
#' calls in a `MultiCompartmapCall` object.
#'
#' @param x A `MultiCompartmapCall` object
#' @param na.omit Whether to omit NAs
#'
#' @concept s7analysis
#' @export
agreement <- new_generic("agreement", "x", function(x, na.omit = TRUE) S7_dispatch())
method(agreement, MultiCompartmapCall) <- function(x, na.omit = TRUE) {
  if (na.omit) {
    get_agreement(na.omit(x@mat))
  } else {
    get_agreement(x@mat)
  }
}

#' Compute correlation between compartment calls
#'
#' Calculates Pearson correlation for every pair of calls in a
#' MultiCompartmapCall object.
#'
#' @param x A `MultiCompartmapCall` object
#' @param na.omit Whether to omit NAs
#'
#' @importFrom stats cor
#' @concept s7analysis
#' @export
corr <- new_generic("corr", "x", function(x, na.omit = TRUE) S7_dispatch())
method(corr, MultiCompartmapCall) <- function(x, na.omit = TRUE) {
  if (na.omit) cor(na.omit(x@mat)) else cor(x@mat)
}

#' Compute the direction of change of singular vectors in the MultiCompartmentCall object
#'
#' The direction of change is computed by chromosome
#'
#' @param x A `MultiCompartmentCall` object
#' @value A `MultiCompartmentCall` object with the gradients in the `@df` slot,
#' replacing the singular vectors
#'
#' @concept s7analysis
#' @export
differentiate <- new_generic("differentiate", "x", function(x) S7_dispatch())
method(differentiate, MultiCompartmapCall) <- function(x) {
  x@df <- x@df[, .(n, pc, name, chr = seqlevels(x))] |>
    _[, .(pc = c(NA, diff(pc[-1])) / diff(n)), by = .(name, chr)] |>
    _[, n := seq_len(.N), by = "name"] |>
    na.omit()
  x@mat <- as.matrix(dcast(x@df, n ~ name, value.var = "pc")[, -1])
  x
}

#' Plot singular values from a `MultiCompartmapCall` object
#'
#' @param x `MultiCompartmapCall`
#' @param ... Placeholder for the `plot` generic - arguments have not effect
#' @param type Whether to plot the singular values as `"line"`plots. Bar plots
#' will be facted by the individual `CompartmapCall` object names while the
#' line plots are overlayed.
#' @param label_coords Label the x-axis with genomic coordinates. Uses a
#' numeric index when set to `FALSE`. Using coordinate labels can severely
#' crowd the x-axis, especially with Kb-resolution calls.
#' @param res The resolution to round the genomic coordinates to (kilobase:
#' "kb" or megabase: "mb")
#' @param width The width of the `geom_line` if `type = "line"` or the width
#' of the bar if `type = "bar"` in the plot
#' @param ylim Upper and lower bound for the y-axis
#'
#' @importFrom ggplot2 ggplot geom_line geom_col scale_y_continuous theme element_text facet_grid vars scale_fill_manual
#' @concept plotting
#' @export
`plot.compartmap::MultiCompartmapCall` <- function(
  x,
  ...,
  type = "line",
  label_coords = FALSE,
  res = "mb",
  width = 1.0,
  ylim = c(-1, 1)
) {
  pd <- x@df
  x_axis <- "n"
  if (label_coords) {
    pd <- x@df[, .(n, pc, name, coord = grscale(x@gr, res))]
    x_axis <- "coord"
  }

  p <- switch(
    type,
    line = {
      ggplot(pd, aes(x = .data[[x_axis]], y = pc, color = name, group = name)) +
        geom_line()
    },
    bar = {
      ggplot(pd, aes(x = .data[[x_axis]], y = pc, group = name, fill = pc > 0)) +
        geom_col(width = width) +
        facet_grid(rows = vars(name)) +
        scale_fill_manual(values = c("deeppink4", "grey50")) +
        theme(legend.position = "none")
    }
  )

  if (label_coords) {
    p <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  p + scale_y_continuous(limits = ylim)
}

#' Get the Mahalanobis distance for the singular values
#' @param mat Matrix of singular values
#' @param cov_method Method to compute covariance. "base": `stats::cov`,
#' "robust": `robust::covRob()`, "mcd": `robust::covRob(estim = "mcd")`
#' @param alpha_level Significance level to use (default: 0.05)
#' @param fdr Whether to perform Benjamini-Hochberg false discovery correction
#'
#' @importFrom stats mahalanobis pchisq
#'
#' @export
compute_mahalanobis <- function(mat, cov_method = c("base", "robust", "mcd"), alpha_level = 0.05, fdr = FALSE) {
  na_row_idx <- which(is.na(mat), arr.ind = TRUE)[, 1]
  norm_mat <- normalize_quantiles(mat)
  full_mat <- na.omit(norm_mat)

  cv <- switch(
    match.arg(cov_method),
    base = cov(full_mat),
    robust = robust::covRob(full_mat)$cov,
    mcd = robust::covRob(full_mat, estim = "mcd")$cov,
  )

  dist <- apply(full_mat, 1, \(i) {
    sqrt(colSums(as.matrix(dist(i)))) / (ncol(full_mat) - 1)
  }) |>
    t()

  dist.mean <- colMeans(dist)
  dist.sd <- apply(dist, 2, sd)
  dist.z <- sweep(dist, 2, dist.mean, FUN = "-") |> sweep(2, dist.sd, FUN = "/")
  pmax <- apply(dist.z, 1, \(i) pnorm(max(i)))
  cen <- full_mat * pmax

  md <- mahalanobis(cen, colMeans(cen), cv)
  pvals <- pchisq(md, df = ncol(mat) - 1, lower.tail = FALSE)
  if (fdr) {
    pvals <- p.adjust(pvals, method = "BH")
  }
  mdf <- data.table(n = seq_len(nrow(mat)))
  mdf[!(n %in% na_row_idx), `:=`(md = md, pval = pvals)][]
}

#' Get bins with significant Mahalanobis distances
#' @param gr Matrix of singular values
#' @param md data.table output from `compute_mahalanobis()`
#' @param alpha_level Significance level to use (default: 0.05)
#' @export
dc_bins <- function(gr, md, alpha_level = 0.05) {
  gr[which(md$pval <= alpha_level)]
}

#' Get the start and end of sequences within a vector of integers
#'
#' It returns a data.table where rows may be duplicated to allow for additional
#' row-specific data if necessary
#'
#' @param v Indices of significant Mahalanobis distances
#' @examples
#' v <- c(1:4, 5, 7:9)
#' get_sequential_idx(v)
#' @export
get_sequential_idx <- function(v) {
  non_consecutive <- v[which(diff(v) != 1)] + 1
  intervals <- findInterval(v, non_consecutive)
  as.data.table(cbind(v, intervals)) |>
    _[, `:=`(start_idx = min(v), end_idx = max(v)), by = intervals] |>
    _[, .(start_idx, end_idx)] |>
    _[, dc_id := .GRP, by = .(start_idx, end_idx)]
}

#' Show differential compartments as overlaid `ggplot2::geom_rect()` over significant distances
#'
#' @param ccall_pd The `@df` slot of a `MultiCompartmapCall` object
#' @param md data.table output from `compute_mahalanobis()`
#' @param alpha_level Significance level to use (default: 0.05)
#' @param show_md Whether to plot the Mahalanobis distance
#' @param fill Color of the `geom_rect`
#' @param alpha Transparency of the `geom_rect`
#' @param ylim The y-axis limits
#' @param xlim The x-axis limits
#' @param label_ids Whether to differential bin ID labels
#' @param select Column name or index to compute differential compartments on.
#' Set this to plot the compartment scores of all input columns but show
#' differential compartments based only on the provided columns.
#'
#' @importFrom ggplot2 ggplot geom_hline geom_rect labs
#' @importFrom patchwork wrap_plots
#' @importFrom stats qchisq
#'
#' @examples
#' v <- c(1:4, 5, 7:9)
#' get_sequential_idx(v)
plot_dc <- function(
  ccall_pd,
  md,
  select,
  type = c("line", "bar"),
  show_md = TRUE,
  alpha_level = 0.05,
  fill = "maroon",
  alpha = 0.5,
  ylim = c(-0.5, 0.5),
  xlim = NULL,
  label_ids = TRUE
) {
  type = match.arg(type)
  pval <- name <- NULL
  xlim <- xlim %||% range(ccall_pd$pos)

  md_pd <- as.data.table(md) |> setnames(c("start", "end"), c("start_idx", "end_idx"))
  dc_pd <- get_sequential_idx(which(md$pval < alpha_level)) |>
    unique() |>
    _[, .(start_pos = start(md[start_idx]), end_pos = end(md[end_idx]), dc_id)]
  cutoff <- qchisq(p = alpha_level, df = length(select) - 1, lower.tail = FALSE)

  cplot <- switch(
    type,
    line = {
      ggplot(ccall_pd, aes(x = pos, y = cscore, color = name)) +
        geom_line() +
        geom_hline(yintercept = 0) +
        scale_y_continuous(limits = ylim) +
        theme(panel.grid = element_blank())
    },
    bar = {
      ggplot(ccall_pd, aes(x = pos, y = cscore, fill = cscore > 0)) +
        geom_col() +
        geom_hline(yintercept = 0) +
        scale_y_continuous(limits = ylim) +
        facet_grid(rows = vars(name)) +
        theme(panel.grid = element_blank())
    }
  )
  chr <- md_pd[, unique(seqnames)]
  cplot <- cplot +
    scale_x_continuous(labels = \(x) x / 1e6, limits = xlim) +
    labs(x = paste(gsub("chr", "Chromosome", chr), "(Mb)"), y = "Compartment score")

  cplot <- cplot +
    geom_rect(
      data = dc_pd,
      stat = "unique",
      inherit.aes = FALSE,
      aes(xmin = start_pos, xmax = end_pos, ymin = ylim[1], ymax = ylim[2]),
      fill = fill,
      alpha = alpha
    )

  if (label_ids) {
    cplot <- cplot +
      geom_label(
        data = dc_pd,
        inherit.aes = FALSE,
        aes(label = dc_id, x = (start_pos + end_pos) / 2, y = ylim[2], vjust = ifelse(dc_id %% 2 == 0, 1.5, 3)),
        size = 2,
        label.size = NA
      )
  }

  if (show_md) {
    mdplot <- as.data.table(md) |>
      setnames(c("start", "end"), c("start_pos", "end_pos")) |>
      ggplot(aes(x = start_pos, y = md)) +
      geom_line() +
      geom_hline(yintercept = cutoff, linetype = "dotted") +
      scale_x_continuous(labels = \(x) x / 1e6, limits = xlim) +
      labs(x = paste(gsub("chr", "Chromosome", chr), "(Mb)"), y = "Mahalanobis distance")

    cplot <- cplot + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    return(wrap_plots(list(cplot, mdplot), nrow = 2))
  }
  cplot
}

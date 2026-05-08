#' Get the Mahalanobis distance for the singular values
#' @param mat Matrix of singular values
#' @param cov_method Method to compute covariance. "base": `stats::cov`,
#' "robust": `robust::covRob()`, "mcd": `robust::covRob(estim = "mcd")`
#'
#' @importFrom stats mahalanobis pchisq
#'
#' @export
compute_mahalanobis <- function(mat, cov_method = c("base", "robust", "mcd")) {
  na_row_idx <- which(is.na(mat), arr.ind = TRUE)[, 1]
  norm_mat <- normalize_quantiles(mat)
  full_mat <- na.omit(norm_mat)

  cv <- switch(
    match.arg(cov_method),
    base = cov(full_mat),
    robust = robust::covRob(full_mat)$cov,
    mcd = robust::covRob(full_mat, estim = "mcd")$cov,
  )

  md <- mahalanobis(full_mat, colMeans(full_mat), cv)
  pvals <- pchisq(md, df = ncol(mat) - 1, lower.tail = FALSE)
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
    _[, `:=`(start = min(v), end = max(v)), by = intervals] |>
    _[, .(idx = v, start, end)]
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
#' @param alpha The significance threshold
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
  type = c("line", "bar"),
  show_md = TRUE,
  alpha_level = 0.05,
  fill = "maroon",
  alpha = 0.5,
  ylim = c(-0.5, 0.5)
) {
  type = match.arg(type)

  pval <- name <- NULL
  seq_idx <- get_sequential_idx(md[, which(pval <= alpha_level)]) |>
    _[, .(start = as.double(start), end = as.double(end))] |>
    unique() |>
    _[start == end, `:=`(start = start - 0.25, end = end + 0.25)]
  cutoff <- qchisq(p = alpha_level, df = ccall_pd[, length(unique(name))] - 1, lower.tail = FALSE)

  cplot <- switch(
    type,
    line = {
      ggplot(ccall_pd, aes(x = n, y = cscore, color = name)) +
        geom_line() +
        geom_hline(yintercept = 0) +
        geom_rect(
          data = seq_idx,
          inherit.aes = FALSE,
          aes(xmin = start, xmax = end, ymin = ylim[1], ymax = ylim[2]),
          fill = fill,
          alpha = alpha
        ) +
        scale_y_continuous(limits = ylim) +
        theme(panel.grid = element_blank())
    },
    bar = {
      ggplot(ccall_pd, aes(x = n, y = cscore, fill = cscore > 0)) +
        geom_col() +
        geom_hline(yintercept = 0) +
        geom_rect(
          data = seq_idx,
          inherit.aes = FALSE,
          aes(xmin = start, xmax = end, ymin = ylim[1], ymax = ylim[2]),
          fill = fill,
          alpha = alpha
        ) +
        scale_y_continuous(limits = ylim) +
        theme(panel.grid = element_blank())
    }
  )

  if (show_md) {
    mdplot <- ggplot(md, aes(x = n, y = md)) +
      geom_line() +
      geom_hline(yintercept = cutoff, linetype = "dotted") +
      labs(y = "Mahalanobis distance")

    return(wrap_plots(list(cplot, mdplot), nrow = 2))
  }
  cplot
}

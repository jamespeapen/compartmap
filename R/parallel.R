#' Check that the number of requested workers is valid
#' @keywords internal
check_worker_count <- function(bpparam, boot.parallel, avail_workers = parallelly::availableCores()) {
  workers <- get_bpnworkers(bpparam)
  total <- sum(Reduce(`*`, workers), workers[1])
  if (verify_workers(total)) {
    return(TRUE)
  }

  msg <- sprintf(
    "Using %1$d outer and %2$d inner workers would require %3$d workers (%1$d + (%1$d x %2$d)) but your system has only %4$d cores.
  See parallelly::availableCores() for more information on available resources",
    workers[1],
    workers[2],
    total,
    avail_workers
  )
  stop(msg)
}

#' Get the total number of BiocParallelParam workers that will get used
#' From a BiocParallelParam or list of 2 BiocParallelParam objects
#' @keywords internal
get_bpnworkers <- function(bp) {
  workers <- bpnworkers.list(bp)
}

#' Return the number of workers in a list of BiocParallelParam objects
#' @param List of BiocParallelParam objects
#' @importFrom BiocParallel bpnworkers
#' @return A vector of the `bpnworkers` count in each list element
#' @keywords internal
bpnworkers.list <- function(bplist) {
  unlist(Map(bpnworkers, bplist))
}

#' Verify that the input BiocParallelParam is valid
#' @param A BiocParallelParam or list of 2 BiocParallelParam objects
#' @importFrom BiocParallel bpnworkers
#' @return TRUE if the total `bpnworkers` in the input does not exceed
#' available resources as defined by `parallelly::availableCores()`
#' @keywords internal
verify_bp <- function(bp) {
  verify_workers(get_bpnworkers(bp))
}

#' Verify that requested thread count is not higher than available
#' @param thread_count The number of workers to check availability
#' @return TRUE if the requested `thread_count` does not exceed available
#' resources as defined by `parallelly::availableCores()`
#' @keywords internal
verify_workers <- function(n_workers) {
  avail_workers <- parallelly::availableCores()
  n_workers <= avail_workers
}

#' Set outer and inner params for nester parallelization
#' The outer param is across the input samples/columns and the second is for
#' bootstrapping. If `boot.parallel` is FALSE, the inner param is set to
#' `SerialParam`.
#' @keywords internal
get_nested_params <- function(BPPARAM, boot.parallel) {
  stopifnot("Only two BiocParallelParam objects can be used" = length(BPPARAM) <= 2)
  single_param <- length(BPPARAM) == 1
  BPPARAM <- if (single_param & is.list(BPPARAM)) BPPARAM[[1]] else BPPARAM

  if (boot.parallel) {
    if (single_param) {
      return(list(outer = BPPARAM, inner = BPPARAM))
    } else {
      return(list(outer = BPPARAM[[1]], inner = BPPARAM[[2]]))
    }
  }

  if (single_param) {
    list(outer = BPPARAM, inner = SerialParam())
  } else {
    list(outer = BPPARAM[[1]], inner = SerialParam())
  }
}

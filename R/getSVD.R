#' Compute the SVD of a matrix using irlba
#'
#' @name getSVD
#'
#' @param matrix A p x n input matrix
#' @param k Number of singular vectors to return
#' @param sing.vec Whether to return the right or left singular vector
#'
#' @return A singular vector or matrix with sign corresponding to positive values
#'
#' @importFrom stats cor
#' @importFrom BiocSingular IrlbaParam runSVD
#' @export
#'
#' @examples
#'
#' dummy <- matrix(rnorm(10000), ncol=25)
#' sing_vec <- getSVD(dummy, k = 1, sing.vec = "right")
getSVD <- function(matrix, k = 1, sing.vec = c("left", "right")) {
  sing.vec <- match.arg(sing.vec)

  matrix <- .centerMatrix(matrix)
  V <- .getUV(matrix, k, sing.vec)
  V <- fixSign(matrix, V)
  V * sqrt(length(V)) # Chromosome length normalization
}

# Centre the matrix subtracting the row means from each row
#' @keywords internal
.centerMatrix <- function(mat) {
  t(scale(t(mat), center = TRUE, scale = FALSE))
}

#' Run SVD and get left or right singular vectors
#'
#' @param matrix A p x n input matrix
#' @param k Number of singular vectors to return
#' @param sing.vec Whether to return the right or left singular vector
#'
#' @return left or right singular vectors
#'
#' @keywords internal
.getUV <- function(matrix, k, sing.vec) {
  SVD <- runSVD(matrix, k = k, BSPARAM = IrlbaParam())
  switch(sing.vec, left = SVD$u, right = SVD$v)
}

#' Fix the sign of SVD singular vectors based on concordance between **V**
#' direction and data direction
#'
#'
#' @details
#'
#' Modified from the reference below. Instead of doing `sign(dot_prod) *
#' (dot_prod^2)` we just take the dot product of the 1st left singular vector
#' **V** and the rows of the correlation matrix. Following the original
#' algorithm exactly results in the values always getting switched to positive
#' which is always what we want here. This also returns a sign 'dominance'
#' metric which the difference between the proportion of positive and negative
#' dot products: a value of 0.8 would mean that 90% of the dot-products were in
#' the dominant direction: 0.9 - 0.1 = 0.8.
#'
#' Resolving the sign ambiguity in the singular value decomposition, Journal of
#' Chemometrics 22(2):135-140
#' (https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/abs/10.1002/cem.1122)
#' Python reference: https://gist.github.com/David-Willo/1825bf9e8c30e13147e332734bcaebd5
#' Discussion of paper and why `abs()` is written in words but not in the algorithm:
#' https://www.mathworks.com/matlabcentral/fileexchange/22118-sign-correction-in-svd-and-pca
#'
#' @param matrix A p x n input matrix
#' @param k Number of singular vectors to return
#' @param sing.vec Whether to return the right or left singular vector
#'
#' @return left or right singular vectors
fixSign <- function(mat, v) {
  dot_prod <- t(v) %*% t(mat) # row-wise
  all_directions <- sign(dot_prod)
  props <- prop.table(table(all_directions))
  direction <- sum(all_directions)
  dominance <- round(abs(props[1] - props[2]), 2) # by how many more were in the dominant direction
  message("Sign dominance: ", dominance)
  if (direction < 0) {
    v <- -v
  }
  v
}

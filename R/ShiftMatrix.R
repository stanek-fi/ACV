#' Construct shift matrix
#'
#' Internal function for creation of sparse shift matrix.
#'
#' @param n Integer specifying dimensions of the shift matrix.
#' @param q Integer specifying the order of the shift matrix. Value `q = 1` (resp. `q = -1`) indicates the upper (resp. lower) shift matrix. Larger (resp. smaller) values represent higher powers of the respective shift matrices.
#'
#' @return Returns a sparse matrix (class `"ngCMatrix"`).
#'
#' @export
#' @keywords internal

shiftMatrix <- function(n, q) {
  Indices <- cbind(
    max(-q, 0) + seq_len(max(n - abs(q), 0)),
    max(q, 0) + seq_len(max(n - abs(q), 0))
  )
  Output <- Matrix::sparseMatrix(i = Indices[, 1], j = Indices[, 2], dims = c(n, n))
  return(Output)
}

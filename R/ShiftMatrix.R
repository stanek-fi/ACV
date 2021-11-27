#' Title
#'
#' Description
#' @export
#' @keywords internal

ShiftMatrix <- function(n, q) {
  Indices <- cbind(
    max(-q, 0) + seq_len(max(n - abs(q), 0)),
    max(q, 0) + seq_len(max(n - abs(q), 0))
  )
  Output <- sparseMatrix(i=Indices[,1],j=Indices[,2],dims=c(n,n))
  return(Output)
}

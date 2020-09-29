#' Title
#'
#' Description
#' @export
#' @keywords internal

ShiftMatrix <- function(n, q) {
  Output <- matrix(0, n, n)
  Indices <- cbind(
    max(-q, 0) + seq_len(max(n - abs(q), 0)),
    max(q, 0) + seq_len(max(n - abs(q), 0))
  )
  Output[Indices] <- 1
  return(Output)
}

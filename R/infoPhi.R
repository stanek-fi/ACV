#' Recover information about `Phi`
#'
#' Recovers all the necessary parameters using which the `Phi` was constructed.
#' @export
#' @keywords internal

infoPhi <- function(Phi) {
  K <- ncol(Phi)
  mn <- nrow(Phi)
  m <- sum(!is.na(Phi[, K]))
  v <- sum(cumprod(is.na(Phi[, 2])))
  h <- sum(!is.na(Phi[, 1])) - m
  mh <- m + h
  J <- lapply(1:K, function(k) 1:sum(!is.na(Phi[, k])))
  return(list(K = K, mn = mn, m = m, v = v, h = h, mh = mh, J = J))
}

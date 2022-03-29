#' Estimate long-run variance
#'
#' Routine for estimating the long-run variance.
#' @export
#' @keywords internal

estimateLongRunVar <- function(x, bw = NULL) {
  if (is.null(bw)) {
    bw <- floor(0.75 * length(x)^(1 / 3))
  }
  weightsACF <- seq(1, 0, by = -(1 / (bw + 1)))[-(bw + 2)]
  ACF <- as.vector(stats::acf(x, type = "covariance", plot = F, lag.max = bw)$acf)
  LongRunVar <- sum(ACF * weightsACF[1:length(ACF)] * c(1, rep(2, bw))[1:length(ACF)])
  return(LongRunVar)
}

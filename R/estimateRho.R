#' Estimate rho coefficient
#'
#' Routine for estimating the rho coefficient.
#' @export
#' @keywords internal

estimateRho <- function(Phi, rhoLimit) {
  # list2env(infoPhi(Phi),environment())
  temp <- infoPhi(Phi)
  K <- temp$K
  mn <- temp$mn
  m <- temp$m
  v <- temp$v
  h <- temp$h
  mh <- temp$mh
  J <- temp$J

  Target <- 1 - 1 / (2 * stats::var(c(Phi),na.rm=T)) * sapply(1:K, function(k) {apply((Phi-Phi[,k])^2,2, mean, na.rm=T)})
  Target <- sapply(1:K, function(k) mean(Target[cbind((k):K, 1:(K - k + 1))]))
  Weights <- pmax(mh - (0:(K - 1)) * v, 0) * (K:1)

  moments <- function(par) {
    return(ifelse(is.na(Target), 0, Target - par^(v * (0:(K - 1)))))
  }

  objective <- function(par) {
    return(sum(moments(par)^2 * Weights))
  }

  fit <- stats::optimize(objective, interval = c(-rhoLimit, rhoLimit))
  rho <- fit$minimum
  return(rho)
}

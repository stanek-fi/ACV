#' Test equality of out-of-sample losses of two algorithms
#'
#' Function `testL()` tests the null hypothesis of equal predictive ability of `algorithm1` and `algorithm2` on time series `y`.
#'
#' @param y Univariate time-series object.
#' @param algorithm1  First algorithm which is to be applied to the time-series. The object which the algorithm produces should respond to fitted and forecast methods.
#' Alternatively in the case of more complex custom algorithms, the algorithm may be a function which takes named arguments `("yInSample", "yOutSample", "h")` or `("yInSample", "yOutSample", "h", "xregInSample", "xregOutSample")` as inputs and produces list with named elements `("yhatInSample", "yhatOutSample")` containing vectors of in-sample and out-of-sample forecasts.
#' @param algorithm2  Second algorithm. See above.
#' @param m Length of the window on which the algorithm should be trained.
#' @param h Number of predictions made after a single training of the algorithm.
#' @param v Number of periods by which the estimation window progresses forward once the predictions are generated.
#' @param xreg Matrix of exogenous regressors supplied to the algorithm (if applicable).
#' @param lossFunction Loss function used to compute contrasts (defaults to squared error).
#' @param method Can be set to either `"optimal"` for the test which optimally utilizes also the in-sample performance or `"convetional"` for the conventional test.
#' @param test Type of the test which is to be executed. Can attain values `"Diebold-Mariano"` for the canonical test of equal predictive ability or `"Ibragimov-Muller"` for the sub-sampling t-test.
#' @param Ha Alternative hypothesis. Can attain values `"!=0"` for two sided test or `"<0"` and `">0"` for one sided tests.
#' @param Phi You can also directly supply `Phi=Phi1-Phi2`; the matrix of contrasts differentials produced by `tsACV`. In this case parameters: `y`, `algorithm`, `m`, `h`, `v`, `xreg`, `lossFunction` are ignored.
#' @param bw Applicable to `"Diebold-Mariano"` test. Bandwidth for the long run variance estimator. If `NULL`, `bw` is selected according to `(3/4)*n^(1/3)`.
#' @param groups  Applicable to `"Ibragimov-Muller"` test. The number of groups to which the data is to be divided.
#' @param ... Other parameters passed to algorithms.
#'
#' @return List containing loss differential estimate and associated p-value along with some other auxiliary information like the matrix of contrasts Phi and the weights used for computation.
#'
#' @examples
#' set.seed(1)
#' y <- rnorm(40)
#' m <- 36
#' h <- 1
#' v <- 1
#' algorithm1 <- function(y) {
#'   forecast::Arima(y, order = c(1, 0, 0))
#' }
#' algorithm2 <- function(y) {
#'   forecast::Arima(y, order = c(2, 0, 0))
#' }
#' testL(y, algorithm1, algorithm2, m = m, h = h, v = v)
#'
#' @export


testL <- function(y, algorithm1, algorithm2, m, h = 1, v = 1, xreg = NULL, lossFunction = function(y, yhat) {(y - yhat)^2}, method = "optimal", test = "Diebold-Mariano", Ha = "!=0", Phi = NULL, bw = NULL, groups = 2, ...) {
  if (is.null(Phi)) {
    Phi1 <- tsACV(y, algorithm1, m, h, v, xreg, lossFunction, ...)
    Phi2 <- tsACV(y, algorithm2, m, h, v, xreg, lossFunction, ...)
    Phi <- Phi1 - Phi2
  }
  temp <- infoPhi(Phi)
  K <- temp$K
  mn <- temp$mn
  m <- temp$m
  v <- temp$v
  h <- temp$h
  mh <- temp$mh
  J <- temp$J
  if (v != h) {
    stop("Currently, only predictive ability testing with h = v is supported.")
  }

  switch(test,
    "Diebold-Mariano" = {
      output <- estimateL(Phi = Phi, method = method, bw = bw)
      estimate <- output$estimate
      tval <- output$estimate / sqrt(output$var)
      pval <- switch(Ha,
        "<0" = stats::pnorm(tval),
        "!=0" = (1 - stats::pnorm(abs(tval))) * 2,
        ">0" = 1 - stats::pnorm(tval)
      )
    },

    "Ibragimov-Muller" = {
      groupSize <- ceiling((K - 1) / groups)
      estimates <- sapply(1:groups, function(g) {
        colIndices <- 1:(groupSize + 1) + (g - 1) * groupSize
        rowIndices <- 1:(m + v * (groupSize)) + (g - 1) * groupSize
        PhiSubset <- Phi[rowIndices[rowIndices <= nrow(Phi)], colIndices[colIndices <= ncol(Phi)]]
        output <- estimateL(Phi = PhiSubset, method = method, bw = bw)
        return(output$estimate)
      })
      estimate <- mean(estimates)
      tval <- estimate / sqrt((1 / (groups)) * (1 / (groups - 1)) * sum((estimates - estimate)^2))
      pval <- switch(Ha,
        "<0" = stats::pt(tval, df = groups - 1),
        "!=0" = (1 - stats::pt(abs(tval), df = groups - 1)) * 2,
        ">0" = 1 - stats::pt(tval, df = groups - 1)
      )
    }
  )

  output <- list(
    estimate = estimate,
    tval = tval,
    pval = pval,
    Ha = Ha,
    test = test,
    Phi = Phi
  )
  class(output) <- "testL"

  return(output)
}

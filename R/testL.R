#' Test equality of losses of two algorithms
#'
#' Function testL test null hypothesis of equal predictive ability of algorithm1 and algorithm2 on time series y. By default, it uses the optimal weighting scheme which exploits also in-sample contrasts in order to deliver more power than regular tests.
#'
#' @param y Univariate time-series object
#' @param algorithm1  First algorithm which is to be applied to the time-series. The object which the algorithm produces should respond to fitted and forecast methods.
#' Alternatively in the case of more complex custom algorithms, the algorithm may be a function which takes named arguments ("yInSample", "yOutSample", "h") or ("yInSample", "yOutSample", "h", "xregInSample", "xregOutSample") as inputs and produces list with named elements ("yhatInSample", "yhatOutSample") containing vectors of in-sample and out-of-sample forecast.
#' @param algorithm2  Second algorithm. See above.
#' @param m Length of the window on which the algorithm ought to be trained.
#' @param h Number of predictions made after single training the algorithm.
#' @param v Number of periods by which the estimation window is shifted once the predictions are generated.
#' @param xreg Matrix of exogenous regressors supplied to the algorithm (if applicable).
#' @param lossFunction Loss function used to compute contrasts (defaults to square loss).
#' @param method Can attain values "augmented"  for the improved estimator which optimally utilizes also in-sample contrast or "regular" for the standard loss estimator.
#' @param test Type of test which is to be executed. Can attain values "Diebold-Mariano" for canonical test of equal predictive ability or "Ibragimov-Muller" sub-sampling t-test.
#' @param Ha Alternative hypothesis. Can attain "!=0" for two sided test or "<0" and ">0" for one sided test.
#' @param Phi One can also directly supply Phi=Phi1-Phi2; the matrix of contrasts differentials produced by tsACV in which case parameters: y, algorithm, m, h, v, xreg, lossFunction are not required.
#' @param bw Applicable to "Diebold-Mariano" test. Bandwidth for long run variance estimator. If null, bw is selected according to common rule of thumb (3/4)*n^(1/3).
#' @param groups  Applicable to "Ibragimov-Muller" test. Number of groups to which the data is to be divided.
#' @param ... Other parameters passed to the algorithm.
#'
#' @return List containing loss differential estimate and associated p-value along with some other auxiliary information like matrix of contrasts Phi and optimal weight vector lambda.
#'
#' @examples
#'
#' mn <- 40
#' y <- rnorm(mn)
#' m <- 36
#' h <- 1
#' v <- 1
#' algorithm1 <- function(y) {
#'   Arima(y, order = c(1, 0, 0))
#' }
#' algorithm2 <- function(y) {
#'   Arima(y, order = c(2, 0, 0))
#' }
#' testL(y, algorithm1, algorithm2, m = m, h = h, v = v)
#'
#' @export


testL <- function(y, algorithm1, algorithm2, m, h = 1, v = 1, xreg = NULL, lossFunction = function(y, yhat) {(y - yhat)^2}, method = "augmented", test = "Diebold-Mariano", Ha = "!=0", Phi = NULL, bw = NULL, groups = 2, ...) {
  if (is.null(Phi)) {
    Phi1 <- tsACV(y, algorithm1, m, h, v, xreg, lossFunction, ...)
    Phi2 <- tsACV(y, algorithm2, m, h, v, xreg, lossFunction, ...)
    Phi <- Phi1 - Phi2
  }
  list2env(infoPhi(Phi), environment())
  if (v != h) {
    stop("Currently, only predictive ability testing with h = v is supported.")
  }

  switch(test,
    "Diebold-Mariano" = {
      output <- estimateL(Phi = Phi, method = method, bw = bw)
      estimate <- output$estimate
      tval <- output$estimate / sqrt(output$var)
      pval <- switch(Ha,
        "<0" = pnorm(tval),
        "!=0" = (1 - pnorm(abs(tval))) * 2,
        ">0" = 1 - pnorm(tval)
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
        "<0" = pt(tval, df = groups - 1),
        "!=0" = (1 - pt(abs(tval), df = groups - 1)) * 2,
        ">0" = 1 - pt(tval, df = groups - 1)
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

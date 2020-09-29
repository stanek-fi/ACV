#' Perform time-series cross-validation for a given algorithm
#'
#' Function tsACV computes contrasts between forecast produced by a given algorithm and the original time-series on which the algorithm is trained.
#' This can be than used to estimate the loss of the algorithm.
#' Unlike the similar tsCV function from forecast package however, taACV also records in-sample contrasts as these can be exploited to produce more accurate loss estimates.
#'
#' @param y Univariate time-series object
#' @param algorithm Algorithm which is to be applied to the time-series. The object which the algorithm produces should respond to fitted and forecast methods.
#' Alternatively in the case of more complex custom algorithms, the algorithm may be a function which takes named arguments ("yInSample", "yOutSample", "h") or ("yInSample", "yOutSample", "h", "xregInSample", "xregOutSample") as inputs and produces list with named elements ("yhatInSample", "yhatOutSample") containing vectors of in-sample and out-of-sample forecast.
#' @param m Length of the window on which the algorithm ought to be trained.
#' @param h Number of predictions made after single training the algorithm.
#' @param v Number of periods by which the estimation window is shifted once the predictions are generated.
#' @param xreg Matrix of exogenous regressors supplied to the algorithm (if applicable).
#' @param lossFunction Loss function used to compute contrasts (defaults to square loss).
#' @param ... Other parameters passed to the algorithm.
#'
#' @return Matrix of computed contrasts. Each row corresponds to particular period of the y time-series and each column corresponds to particular location of the training window.
#'
#' @examples
#'
#' mn <- 40
#' y <- rnorm(mn)
#' m <- 36
#' h <- 1
#' v <- 1
#' algorithm <- function(y) {
#'   Arima(y, order = c(1, 0, 0))
#' }
#' Phi <- tsACV(y, algorithm, m = m, h = h, v = v)
#' Phi
#'
#' @export

tsACV <- function(y, algorithm, m, h = 1, v = 1, xreg = NULL, lossFunction = function(y, yhat) {(y - yhat)^2}, ...) {
  mn <- length(y)
  y <- ts(c(y, y[rep(1, h)]), start = start(y), frequency = frequency(y))
  if (!is.null(xreg)) {
    if (!ifelse(is.matrix(xreg), (nrow(xreg) == mn), F)) {
      stop("xreg must be a matrix with the same number of rows as length(y)")
    }
    xreg <- ts(rbind(xreg, xreg[rep(1, h), ]), start = start(xreg), frequency = frequency(xreg))
  }
  if (((mn - m) <= 0) | ((mn - m) %% v != 0)) {
    stop("number of periods left for evaluation must be greater than 0 and divisble by h")
  }
  I <- seq(0, (mn - m), by = v)
  yhat <- ts(matrix(NA_real_, nrow = mn, ncol = length(I)))

  for (index in seq_along(I)) {
    i <- I[index]
    yInSample <- window(y, 1 + i, m + i)
    yOutSample <- window(y, 1 + m + i, h + m + i)

    if (is.null(xreg)) {
      if(all(c("yInSample", "yOutSample", "h")%in%formalArgs(algorithm))){
        model <- algorithm(yInSample = yInSample, yOutSample = yOutSample, h = h, ...)
        yhatInSample <- model$yhatInSample
        yhatOutSample <- model$yhatOutSample
      }else{
        model <- algorithm(yInSample, ...)
        yhatOutSample <- forecast(model, h = h)$mean
        yhatInSample <- fitted(model)
      }
    } else {
      xregInSample <- window(xreg, 1 + i, m + i)
      xregOutSample <- window(xreg, 1 + m + i, h + m + i)

      if(all(c("yInSample", "yOutSample", "h", "xregInSample", "xregOutSample")%in%formalArgs(algorithm))){
        model <- algorithm(yInSample = yInSample, yOutSample = yOutSample, h = h, xregInSample = xregInSample, xregOutSample = xregOutSample, ...)
        yhatInSample <- model$yhatInSample
        yhatOutSample <- model$yhatOutSample
      }else{
        model <- algorithm(yInSample, xreg = xregInSample, ...)
        yhatOutSample <- forecast(model, xreg = xregOutSample, h = h)$mean
        yhatInSample <- fitted(model)
      }
    }

    timeindices <- (1:(m + h) + i)
    yhat[timeindices[timeindices<=mn], index] <- c(yhatInSample, yhatOutSample)[timeindices<=mn]
  }
  output <- apply(yhat, 2, function(x) lossFunction(y[seq_len(mn)], x))
  colnames(output) <- paste("i=", I, sep = "")
  return(output)
}

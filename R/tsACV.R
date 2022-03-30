#' Perform time-series cross-validation
#'
#' Function `tsACV()` computes contrasts between forecasts produced by a given algorithm and the original time-series on which the algorithm is trained.
#' This can then be used to estimate the loss of the algorithm.
#' Unlike the similar `tsCV()` function from the `'forecast'` package, `tsACV()` also records in-sample contrasts as these can be leveraged to produce more accurate out-of-sample loss estimates.
#'
#' @param y Univariate time-series object.
#' @param algorithm Algorithm which is to be applied to the time-series. The object which the algorithm produces should respond to fitted and forecast methods.
#' Alternatively in the case of more complex custom algorithms, the algorithm may be a function which takes named arguments `("yInSample", "yOutSample", "h")` or `("yInSample", "yOutSample", "h", "xregInSample", "xregOutSample")` as inputs and produces list with named elements `("yhatInSample", "yhatOutSample")` containing vectors of in-sample and out-of-sample forecasts.
#' @param m Length of the window on which the algorithm should be trained.
#' @param h Number of predictions made after a single training of the algorithm.
#' @param v Number of periods by which the estimation window progresses forward once the predictions are generated.
#' @param xreg Matrix of exogenous regressors supplied to the algorithm (if applicable).
#' @param lossFunction Loss function used to compute contrasts (defaults to squared error).
#' @param ... Other parameters passed to the algorithm.
#'
#' @return Matrix of computed contrasts. Each row corresponds to a particular period of the y time-series and each column corresponds to a particular location of the training window.
#'
#' @examples
#' set.seed(1)
#' y <- rnorm(40)
#' m <- 36
#' h <- 1
#' v <- 1
#' tsACV(y, forecast::Arima, m = m, h = h, v = v)
#'
#' @export

tsACV <- function(y, algorithm, m, h = 1, v = 1, xreg = NULL, lossFunction = function(y, yhat) {(y - yhat)^2}, ...) {
  mn <- length(y)
  if(any(is.na(y))){
    stop("y should not contain missing values")
  }
  y <- stats::ts(c(y, y[rep(1, h)]), start = stats::start(y), frequency = stats::frequency(y))
  if (!is.null(xreg)) {
    if (!ifelse(is.matrix(xreg), (nrow(xreg) == mn), F)) {
      stop("xreg must be a matrix with the same number of rows as length(y)")
    }
    xreg <- stats::ts(rbind(xreg, xreg[rep(1, h), ]), start = stats::start(xreg), frequency = stats::frequency(xreg))
  }
  if (((mn - m) <= 0) | ((mn - m) %% v != 0)) {
    stop("number of periods left for evaluation must be greater than 0 and divisble by h")
  }
  I <- seq(0, (mn - m), by = v)
  yhat <- stats::ts(matrix(NA_real_, nrow = mn, ncol = length(I)))
  times <- stats::time(y)

  for (index in seq_along(I)) {
    i <- I[index]
    yInSample <- stats::window(y, times[1 + i], times[m + i])
    yOutSample <- stats::window(y, times[1 + m + i], times[h + m + i])

    if (is.null(xreg)) {
      if(all(c("yInSample", "yOutSample", "h")%in%methods::formalArgs(algorithm))){
        model <- algorithm(yInSample = yInSample, yOutSample = yOutSample, h = h, ...)
        yhatInSample <- model$yhatInSample
        yhatOutSample <- model$yhatOutSample
      }else{
        model <- algorithm(yInSample, ...)
        yhatOutSample <- forecast::forecast(model, h = h)$mean
        yhatInSample <- stats::fitted(model)
      }
    } else {
      xregInSample <- stats::window(xreg, times[1 + i], times[m + i])
      xregOutSample <- stats::window(xreg, times[1 + m + i], times[h + m + i])

      if(all(c("yInSample", "yOutSample", "h", "xregInSample", "xregOutSample")%in%methods::formalArgs(algorithm))){
        model <- algorithm(yInSample = yInSample, yOutSample = yOutSample, h = h, xregInSample = xregInSample, xregOutSample = xregOutSample, ...)
        yhatInSample <- model$yhatInSample
        yhatOutSample <- model$yhatOutSample
      }else{
        model <- algorithm(yInSample, xreg = xregInSample, ...)
        yhatOutSample <- forecast::forecast(model, xreg = xregOutSample, h = h)$mean
        yhatInSample <- stats::fitted(model)
      }
    }

    timeindices <- (1:(m + h) + i)
    yhat[timeindices[timeindices<=mn], index] <- c(yhatInSample, yhatOutSample)[timeindices<=mn]
  }
  output <- apply(yhat, 2, function(x) lossFunction(y[seq_len(mn)], x))
  colnames(output) <- paste("i=", I, sep = "")
  return(output)
}

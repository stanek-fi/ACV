#' Estimate loss of a given algorithm
#'
#' Function estimateL estimate loss of a given algorithm on time-series y. By default, it uses the optimal weighting scheme which exploits also in-sample contrasts in order to deliver more precise estimates than regular estimators.
#'
#' @param y Univariate time-series object
#' @param algorithm Algorithm which is to be applied to the time-series. The object which the algorithm produces should respond to fitted and forecast methods.
#' Alternatively in the case of more complex custom algorithms, the algorithm may be a function which takes named arguments ("yInSample", "yOutSample", "h") or ("yInSample", "yOutSample", "h", "xregInSample", "xregOutSample") as inputs and produces list with named elements ("yhatInSample", "yhatOutSample") containing vectors of in-sample and out-of-sample forecast.
#' @param m Length of the window on which the algorithm ought to be trained.
#' @param h Number of predictions made after single training the algorithm.
#' @param v Number of periods by which the estimation window is shifted once the predictions are generated.
#' @param xreg Matrix of exogenous regressors supplied to the algorithm (if applicable).
#' @param lossFunction Loss function used to compute contrasts (defaults to square loss).
#' @param method Can attain values "augmented"  for the improved estimator which optimally utilizes also in-sample contrast or "regular" for the standard loss estimator.
#' @param Phi One can also directly supply Phi; the matrix of contrasts produced by tsACV in which case parameters: y, algorithm, m, h, v, xreg, lossFunction are not required.
#' @param bw Bandwidth for long run variance estimator. If null, bw is selected according to common rule of thumb (3/4)*n^(1/3).
#' @param ... Other parameters passed to the algorithm.
#'
#' @return List containing loss estimate and its estimated variance along with some other auxiliary information like matrix of contrasts Phi and optimal weight vector lambda.
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
#' estimateL(y, algorithm, m = m, h = h, v = v)
#'
#' @export

estimateL <- function(y, algorithm, m, h = 1, v = 1, xreg = NULL, lossFunction = function(y, yhat) {(y - yhat)^2}, method = "augmented", Phi = NULL, bw = NULL, rhoLimit = 1, ...) {

  if (is.null(Phi)) {
    Phi <- tsACV(y, algorithm, m, h, v, xreg, lossFunction, ...)
  }
  list2env(infoPhi(Phi), environment())

  switch(method,

         "regular" = {
           lambda <- do.call(c, J) > m
           lambda <- lambda / sum(lambda)
           if (v != h) {
             warning("Currently, standart error estimation is supported only for h = v.")
             var <- NA
           } else {
             phiOutSample <- na.omit(c(Phi))[do.call(c, J) > m]
             var <- estimateLongRunVar(phiOutSample, bw) / length(phiOutSample)
           }
           rho <- NA
         },

         "augmented" = {
           b <- sapply(1:mh, function(x) {
             ifelse(x > m, sum(sapply(J, function(Jk) {
               x %in% Jk
             })), 0)
           })
           b <- b / sum(b)
           rho <- estimateRho(Phi, rhoLimit)

           I <- ShiftMatrix(mh, 0)
           Z1 <- ShiftMatrix(mh, 0) + rho^2 / (1 - rho^2) * ShiftMatrix(mh, -v) %*% ShiftMatrix(mh, v)
           Z2 <- ShiftMatrix(mh, 0) + rho^2 / (1 - rho^2) * (ShiftMatrix(mh, -v) %*% ShiftMatrix(mh, v) + ShiftMatrix(mh, v) %*% ShiftMatrix(mh, -v))
           Z3 <- ShiftMatrix(mh, 0) + rho^2 / (1 - rho^2) * ShiftMatrix(mh, v) %*% ShiftMatrix(mh, -v)
           Zu <- -rho / (1 - rho^2) * ShiftMatrix(mh, v)
           Zl <- -rho / (1 - rho^2) * ShiftMatrix(mh, -v)

           BViB <- I[, J[[1]]] %*% Z1[J[[1]], J[[1]]] %*% I[J[[1]], ] + I[, J[[2]]] %*% Zu[J[[2]], J[[1]]] %*% I[J[[1]], ]
           for (k in seq(2, K - 1, length = max(0, K - 2))) {
             BViB <- BViB + I[, J[[k - 1]]] %*% Zl[J[[k - 1]], J[[k]]] %*% I[J[[k]], ] + I[, J[[k]]] %*% Z2[J[[k]], J[[k]]] %*% I[J[[k]], ] + I[, J[[k + 1]]] %*% Zu[J[[k + 1]], J[[k]]] %*% I[J[[k]], ]
           }
           BViB <- BViB + I[, J[[K - 1]]] %*% Zl[J[[K - 1]], J[[K]]] %*% I[J[[K]], ] + I[, J[[K]]] %*% Z3[J[[K]], J[[K]]] %*% I[J[[K]], ]
           BViBi <- solve(BViB)

           ViB <- do.call(rbind, c(
             list(Z1[J[[1]], J[[1]]] %*% I[J[[1]], ] + Zl[J[[1]], J[[2]]] %*% I[J[[2]], ]),
             lapply(seq(2, K - 1, length = max(0, K - 2)), function(k) {
               Zu[J[[k]], J[[k - 1]]] %*% I[J[[k - 1]], ] + Z2[J[[k]], J[[k]]] %*% I[J[[k]], ] + Zl[J[[k]], J[[k + 1]]] %*% I[J[[k + 1]], ]
             }),
             list(Zu[J[[K]], J[[K - 1]]] %*% I[J[[K - 1]], ] + Z3[J[[K]], J[[K]]] %*% I[J[[K]], ])
           ))

           lambda <- as.vector(ViB %*% BViBi %*% b)

           if(v != h){
             warning("Currently, standart error estimation is supported only for h = v.")
             var <- NA
           }else{
             phiOutSample <- na.omit(c(Phi))[do.call(c, J) > m]
             varRatio <- c(t(b) %*% BViBi %*% b) / (1 / length(phiOutSample))
             var <- estimateLongRunVar(phiOutSample, bw) / length(phiOutSample) * varRatio
           }
         }
  )

  output <- list(
    estimate = sum(na.omit(c(Phi)) * lambda),
    var = var,
    lambda = lambda,
    Phi = Phi,
    rho = rho
  )
  class(output) <- "estimateL"
  return(output)
}

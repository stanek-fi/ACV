% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testL.R
\name{testL}
\alias{testL}
\title{Test equality of losses of two algorithms}
\usage{
testL(
  y,
  algorithm1,
  algorithm2,
  m,
  h = 1,
  v = 1,
  xreg = NULL,
  lossFunction = function(y, yhat) {     (y - yhat)^2 },
  method = "augmented",
  test = "Diebold-Mariano",
  Ha = "!=0",
  Phi = NULL,
  bw = NULL,
  groups = 2,
  ...
)
}
\arguments{
\item{y}{Univariate time-series object}

\item{algorithm1}{First algorithm which is to be applied to the time-series. The object which the algorithm produces should respond to fitted and forecast methods.
Alternatively in the case of more complex custom algorithms, the algorithm may be a function which takes named arguments ("yInSample", "yOutSample", "h") or ("yInSample", "yOutSample", "h", "xregInSample", "xregOutSample") as inputs and produces list with named elements ("yhatInSample", "yhatOutSample") containing vectors of in-sample and out-of-sample forecast.}

\item{algorithm2}{Second algorithm. See above.}

\item{m}{Length of the window on which the algorithm ought to be trained.}

\item{h}{Number of predictions made after single training the algorithm.}

\item{v}{Number of periods by which the estimation window is shifted once the predictions are generated.}

\item{xreg}{Matrix of exogenous regressors supplied to the algorithm (if applicable).}

\item{lossFunction}{Loss function used to compute contrasts (defaults to square loss).}

\item{method}{Can attain values "augmented"  for the improved estimator which optimally utilizes also in-sample contrast or "regular" for the standard loss estimator.}

\item{test}{Type of test which is to be executed. Can attain values "Diebold-Mariano" for canonical test of equal predictive ability or "Ibragimov-Muller" sub-sampling t-test.}

\item{Ha}{Alternative hypothesis. Can attain "!=0" for two sided test or "<0" and ">0" for one sided test.}

\item{Phi}{One can also directly supply Phi=Phi1-Phi2; the matrix of contrasts differentials produced by tsACV in which case parameters: y, algorithm, m, h, v, xreg, lossFunction are not required.}

\item{bw}{Applicable to "Diebold-Mariano" test. Bandwidth for long run variance estimator. If null, bw is selected according to common rule of thumb (3/4)*n^(1/3).}

\item{groups}{Applicable to "Ibragimov-Muller" test. Number of groups to which the data is to be divided.}

\item{...}{Other parameters passed to the algorithm.}
}
\value{
List containing loss differential estimate and associated p-value along with some other auxiliary information like matrix of contrasts Phi and optimal weight vector lambda.
}
\description{
Function testL test null hypothesis of equal predictive ability of algorithm1 and algorithm2 on time series y. By default, it uses the optimal weighting scheme which exploits also in-sample contrasts in order to deliver more power than regular tests.
}
\examples{

mn <- 40
y <- rnorm(mn)
m <- 36
h <- 1
v <- 1
algorithm1 <- function(y) {
  Arima(y, order = c(1, 0, 0))
}
algorithm2 <- function(y) {
  Arima(y, order = c(2, 0, 0))
}
testL(y, algorithm1, algorithm2, m = m, h = h, v = v)

}

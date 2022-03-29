library(testthat)
library(ACV)

set.seed(1)
mn <- 20
y <- rnorm(mn)
xreg <- cbind(y + rnorm(mn, 0, 0.1), rnorm(mn))
m <- 14
h <- 2
v <- 1

algorithm <- function(yInSample, yOutSample, h, xregInSample, xregOutSample) {
  dataInSample <- data.frame(y = yInSample, xregInSample)
  dataOutSample <- data.frame(y = NA, xregOutSample)
  model <- lm(y ~ ., dataInSample)
  return(list(
    yhatInSample = fitted(model),
    yhatOutSample = predict(model, dataOutSample)
  ))
}

test_that("tsACV", {
  Phi1 <- tsACV(y, algorithm, m, h = h, v = v, xreg = xreg)
  algorithm <- forecast::Arima
  Phi2 <- tsACV(y, algorithm, m, h = h, v = v, xreg = xreg)
  expect_true(mean(abs(Phi1 - Phi2), na.rm = T) < 1e-10)
})

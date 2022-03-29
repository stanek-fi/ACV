library(testthat)
library(ACV)

mn <- 40
start=c(2000,1)
frequency=4
algorithm <- forecast::Arima
m <- 30
h <- 1
v <- 1

test_that("tsACV under ts input with frequency", {
  set.seed(1)
  y <- rnorm(mn)
  y=ts(y,start=start, frequency = frequency)
  xreg <- cbind(y + rnorm(mn, 0, 0.1), rnorm(mn))
  Phi1 <- tsACV(y, algorithm, m, h = h, v = v, xreg = xreg)

  set.seed(1)
  y <- rnorm(mn)
  xreg <- cbind(y + rnorm(mn, 0, 0.1), rnorm(mn))
  Phi2 <- tsACV(y, algorithm, m, h = h, v = v, xreg = xreg)

  expect_true(mean(abs(Phi1 - Phi2), na.rm = T) < 1e-10)
})






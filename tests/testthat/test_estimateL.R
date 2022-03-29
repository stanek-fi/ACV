library(testthat)
library(ACV)

set.seed(1)
mn <- 40
y <- rnorm(mn)
xreg <- cbind(y + rnorm(mn, 0, 0.1), rnorm(mn))
algorithm <- forecast::Arima
m <- 35
h <- 1
v <- 1
# h <- 5
# v <- 5
# index=1
Phi <- tsACV(y, algorithm, m, h = h, v = v, xreg = xreg)

test_that("estimateLACV", {
  LACVhat1 <- estimateL(y, algorithm, m, h, v, xreg, method = "optimal", rhoLimit = 1)
  LACVhat2 <- estimateL(Phi = Phi, method = "optimal", rhoLimit = 1)
  expect_true(abs(LACVhat1$estimate - 0.006953053) < 1e-9)
  expect_true(abs(LACVhat1$estimate - LACVhat2$estimate) < 1e-10)
})

test_that("estimateLCV", {
  LCVhat1 <- estimateL(y, algorithm, m, h, v, xreg, method = "conventional")
  LCVhat2 <- estimateL(Phi = Phi, method = "conventional")
  expect_true(abs(LCVhat1$estimate - LCVhat2$estimate) < 1e-10)
})

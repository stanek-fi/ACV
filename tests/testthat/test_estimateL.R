library(testthat)
library(ACV)

set.seed(1)
mn <- 40
y <- rnorm(mn)
xreg <- cbind(y + rnorm(mn, 0, 0.1), rnorm(mn))
algorithm <- Arima
m <- 35
h <- 1
v <- 1
Phi <- tsACV(y, algorithm, m, h = h, v = v, xreg = xreg)

test_that("estimateLACV", {
  LACVhat1 <- estimateL(y, algorithm, m, h, v, xreg, method = "augmented")
  LACVhat2 <- estimateL(Phi = Phi, method = "augmented")
  expect_true(abs(LACVhat1$estimate - LACVhat2$estimate) < 1e-10)
})

test_that("estimateLCV", {
  LCVhat1 <- estimateL(y, algorithm, m, h, v, xreg, method = "regular")
  LCVhat2 <- estimateL(Phi = Phi, method = "regular")
  expect_true(abs(LCVhat1$estimate - LCVhat2$estimate) < 1e-10)
})

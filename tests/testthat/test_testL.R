library(testthat)
library(ACV)


set.seed(1)
mn <- 40
y <- rnorm(mn)
m <- 30
h <- 1
v <- 1
algorithm1 <- function(y) {
  forecast::Arima(y, order = c(1, 0, 0))
}
algorithm2 <- function(y) {
  forecast::Arima(y, order = c(2, 0, 0))
}


test_that("testL", {
  regularEst1 <- estimateL(y, algorithm1, m = m, h = h, v = v, method = "conventional")
  augmentedEst1 <- estimateL(y, algorithm1, m = m, h = h, v = v, method = "optimal")
  Phi1 <- tsACV(y, algorithm1, m = m, h = h, v = v)
  regularEst1Alt <- estimateL(Phi = Phi1, method = "conventional")
  augmentedEst1Alt <- estimateL(Phi = Phi1, method = "optimal")

  expect_true(identical(regularEst1$estimate, regularEst1Alt$estimate))
  expect_true(identical(augmentedEst1$estimate, augmentedEst1Alt$estimate))

  regularEst2 <- estimateL(y, algorithm2, m = m, h = h, v = v, method = "conventional")
  augmentedEst2 <- estimateL(y, algorithm2, m = m, h = h, v = v, method = "optimal")
  Phi2 <- tsACV(y, algorithm2, m = m, h = h, v = v)
  regularEst2Alt <- estimateL(Phi = Phi2, method = "conventional")
  augmentedEst2Alt <- estimateL(Phi = Phi2, method = "optimal")

  expect_true(identical(regularEst2$estimate, regularEst2Alt$estimate))
  expect_true(identical(augmentedEst2$estimate, augmentedEst2Alt$estimate))

  regularDMTest12 <- testL(y, algorithm1, algorithm2, m = m, h = h, v = v, method = "conventional", test = "Diebold-Mariano")
  regularDMTest12Alt <- testL(Phi = Phi1 - Phi2, method = "conventional", test = "Diebold-Mariano")
  regularIMTest12 <- testL(y, algorithm1, algorithm2, m = m, h = h, v = v, method = "conventional", test = "Ibragimov-Muller")
  regularIMTest12Alt <- testL(Phi = Phi1 - Phi2, method = "conventional", test = "Ibragimov-Muller")
  augmentedDMTest12 <- testL(y, algorithm1, algorithm2, m = m, h = h, v = v, method = "optimal", test = "Diebold-Mariano")
  augmentedDMTest12Alt <- testL(Phi = Phi1 - Phi2, method = "optimal", test = "Diebold-Mariano")
  augmentedIMTest12 <- testL(y, algorithm1, algorithm2, m = m, h = h, v = v, method = "optimal", test = "Ibragimov-Muller")
  augmentedIMTest12Alt <- testL(Phi = Phi1 - Phi2, method = "optimal", test = "Ibragimov-Muller")

  expect_true(identical(regularDMTest12, regularDMTest12Alt))
  expect_true(identical(regularIMTest12, regularIMTest12Alt))
  expect_true(identical(augmentedDMTest12, augmentedDMTest12Alt))
  expect_true(identical(augmentedIMTest12, augmentedIMTest12Alt))
  expect_true((regularDMTest12$estimate - (regularEst1$estimate - regularEst2$estimate)) < 1e-10)
  expect_true((regularIMTest12$estimate - (regularEst1$estimate - regularEst2$estimate)) < 1e-10)
})

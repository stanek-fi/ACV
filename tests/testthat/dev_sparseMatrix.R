library(testthat)
library(ACV)


set.seed(1)
mn <- 500
mn <- 1000
y <- rnorm(mn)
m <- 450
m <- 950
h <- 1
v <- 1
algorithm1 <- function(y) {
  Arima(y, order = c(1, 0, 0))
}


Phi <- tsACV(y, algorithm1, m = m, h = h, v = v)
start <- Sys.time()
res <- estimateL(Phi = Phi, method = "augmented")
Sys.time() - start



library(Matrix)
library(SparseM)




ShiftMatrix <- function(n, q) {
  # Output <- matrix(0, n, n)
  Indices <- cbind(
    max(-q, 0) + seq_len(max(n - abs(q), 0)),
    max(q, 0) + seq_len(max(n - abs(q), 0))
  )
  # Output[Indices] <- 1
  Output <- sparseMatrix(i=Indices[,1],j=Indices[,2],dims=c(n,n))
  return(Output)
}



list2env(infoPhi(Phi), environment())


start <- Sys.time()
for(interation in 1:10){
  b <- sapply(1:mh, function(x) {
    ifelse(x > m, sum(sapply(J, function(Jk) {
      x %in% Jk
    })), 0)
  })
  b <- b / sum(b)
  rho <- estimateRho(Phi)

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

  sum(na.omit(c(Phi)) * lambda) - res$estimate
}
Sys.time() - start


















# n=1000
# q=1
#
#
# Output <- matrix(0, n, n)
# Indices <- cbind(
#   max(-q, 0) + seq_len(max(n - abs(q), 0)),
#   max(q, 0) + seq_len(max(n - abs(q), 0))
# )
# Output[Indices] <- 1
#
# vec=as.matrix(1:n)
#
#
# M <- Output
# start <- Sys.time()
# t(vec) %*% (3 * M) %*% M %*% vec
# Sys.time() - start
#
#
# M <- sparseMatrix(i=Indices[,1],j=Indices[,2],dims=c(n,n))
# start <- Sys.time()
# t(vec) %*% (3 * M) %*% M %*% vec
# Sys.time() - start
#
#
# M <- as.matrix.csr(Output)
# start <- Sys.time()
# t(vec) %*% (3 * M) %*% M %*% vec
# Sys.time() - start







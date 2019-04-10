plot(1:1000, a$S[, 2], type = "l", xlab = "", ylab = "")
#-----------------------------------------------------------
#Example 3: using FastICA to perform projection pursuit on a
# mixture of bivariate normal distributions
#-----------------------------------------------------------
if(require(MASS)){
  x <- mvrnorm(n = 1000, mu = c(0, 0), Sigma = matrix(c(10, 3, 3, 1), 2, 2))
  x1 <- mvrnorm(n = 1000, mu = c(-1, 2), Sigma = matrix(c(10, 3, 3, 1), 2, 2))
  X <- rbind(x, x1)
  a <- fastICA(X, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
               method = "R", row.norm = FALSE, maxit = 200,
               tol = 0.0001, verbose = TRUE)
  par(mfrow = c(1, 3))
  plot(a$X, main = "Pre-processed data")
  plot(a$X %*% a$K, main = "PCA components")
  plot(a$S, main = "ICA components")
}


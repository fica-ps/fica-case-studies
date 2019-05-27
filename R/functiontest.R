library(fastICA)

# Source matrix
S <- cbind(sin((1:1000)/20), rep((((1:200)-100)/100), 5))
# Mixing matrix
A <- matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
# plot graphs
par(mfcol = c(1, 2))
plot(1:1000, S[,1], type = "l",xlab = "S1", ylab = "")
plot(1:1000, S[,2], type = "l", xlab = "S2", ylab = "")


# Mixed two signals
X <- S %*% A

par(mfcol = c(1, 2))
plot(1:1000, X[,1], type = "l",xlab = "X1", ylab = "")
plot(1:1000, X[,2], type = "l", xlab = "X2", ylab = "")

#initw <- matrix(c(1.02081, 0.408655, -1.92523 ,-0.756068), 2, 2 , TRUE)
initw<- matrix(rnorm(2^2),2,2)
#PCA Whitening to match Julia implementation
S = cov(X)
W = whiteningMatrix (S, method="PCA")
X1 = tcrossprod(X, W) # whitened data


# ICA for extracting independent sources from mixed signals
a <- ica.R.def(t(X1), 2, fun = "logcosh", alpha = 1, maxit = 200,
             tol = 0.0001, verbose = TRUE, initw)
K = a%*%W
S = K%*%t(X)
S = t(S)
par(mfcol = c(1, 2))
plot(1:1000, S[,1], type = "l", xlab = "S'1", ylab = "")
plot(1:1000, S[,2], type = "l", xlab = "S'2", ylab = "")
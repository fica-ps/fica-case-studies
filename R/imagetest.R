
library(imager)



#cat_one
setwd(dirname(sys.frame(1)$ofile))

c_one <- load.image('C:/Users/David/Desktop/ISEL/Projecto/fica-case-studies/R/cat_one_resize.jpg')

c_one.g <- grayscale(c_one)
plot(c_one.g) #Cat one

#cat_two

c_two <- load.image('C:/Users/David/Desktop/ISEL/Projecto/fica-case-studies/R/cat_two_resize.jpg')

c_two.g <- grayscale(c_two) 
plot(c_two.g) #Cat two

#mix both images

rows = nrow(c_one.g)
cols = ncol(c_one.g)

vec_one <- as.vector(array(c_one.g))
vec_two <- as.vector(array(c_two.g))

X = vec_one * 0.6 + vec_two * 0.4

plot(as.cimg(matrix(X,rows,cols)))

Z = vec_one * 0.4 + vec_two * 0.6

plot(as.cimg(matrix(Z,rows,cols)))

X <- cbind(X,Z)

initw <- matrix(c(0.86456941, 1.0646533 , 0.07053282 ,0.8470662), 2, 2 , TRUE)
#initw<- matrix(rnorm(2^2),2,2)
#PCA Whitening to match Julia implementation
S = cov(X)
W = whiteningMatrix (S, method="PCA")
X1 = tcrossprod(X, W) # whitened data

a <- ica.R.def(t(X1), 2, fun = "logcosh", alpha = 1, maxit = 200,
               tol = 0.0001, verbose = TRUE, initw)
K = a%*%W
S = X%*%t(K)
S = t(S)


#a <- fastICA(K, 2, alg.typ = "deflation", fun = "logcosh", alpha = 1,
             #method = "R", row.norm = FALSE, maxit = 100,
             #tol = 0.0001, verbose = TRUE)
plot(as.cimg(matrix(S[1,]+S[2,],rows,cols)), main = "Cat one")
plot(as.cimg(matrix(S[1,]-S[2,],rows,cols)), main = "Cat two")
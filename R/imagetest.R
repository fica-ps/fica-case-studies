
library(imager)



#cat_one
setwd(dirname(sys.frame(1)$ofile))

c_one <- load.image('C:/Users/leonw/Desktop/ISEL/project/cran/fastICA/R/cat_one_resize.jpg')

c_one.g <- grayscale(c_one)
plot(c_one.g) #Cat one

#cat_two

c_two <- load.image('C:/Users/leonw/Desktop/ISEL/project/cran/fastICA/R/cat_two_resize.jpg')

c_two.g <- grayscale(c_two) 
plot(c_two.g) #Cat two

#mix both images

rows = nrow(c_one.g)
cols = ncol(c_one.g)

vec_one <- as.vector(array(c_one.g))
vec_two <- as.vector(array(c_two.g))

X = vec_one * 0.2 + vec_two * 0.8

plot(as.cimg(matrix(X,rows,cols)))

Z = vec_one * 0.8 + vec_two * 0.2

plot(as.cimg(matrix(Z,rows,cols)))

K <- rbind(X,Z)

a <- fastICA(K, 2, alg.typ = "deflation", fun = "logcosh", alpha = 1,
             method = "R", row.norm = FALSE, maxit = 100,
             tol = 0.0001, verbose = TRUE)

plot(as.cimg(matrix(a$S[,1],rows,cols)), main = "Cat one")
plot(as.cimg(matrix(a$S[,2],rows,cols)), main = "Cat two")
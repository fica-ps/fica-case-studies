### whiteningMatrix.R  (2019-01-09)
###
###    Compute whitening matrix
###
### Copyright 2018-19 Korbinian Strimmer
###
###
### This file is part of the `whitening' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


# create whitening matrix W from given covariance matrix Sigma
whiteningMatrix = function(Sigma, 
  method=c("ZCA", "PCA", "Cholesky", "ZCA-cor", "PCA-cor"))
{
  method=match.arg(method)
  if(method=="ZCA" | method=="PCA")
  {
    eSigma = eigen(Sigma, symmetric=TRUE) 
    U = eSigma$vectors
    lambda = eSigma$values
    
    # fix sign ambiguity in eigenvectors by making U positive diagonal
    #U = sweep(U, 2, sign(diag(U)), "*") # U %*% diag( sign(diag(U)) )

    W = diag(1/sqrt(lambda)) %*% t(U)
    if (method=="ZCA") W = U %*% W
  }

  if(method=="Cholesky")
  {
     W = chol(solve(Sigma))
  }

  if(method=="ZCA-cor" | method=="PCA-cor")
  {
    v = diag(Sigma)
    R = cov2cor(Sigma)
    eR = eigen(R, symmetric=TRUE) 
    G = eR$vectors
    theta = eR$values

    # fix sign ambiguity in eigenvectors by making G positive diagonal
    G = sweep(G, 2, sign(diag(G)), "*") # G %*% diag( sign(diag(G)) )

    W = diag(1/sqrt(theta)) %*% t(G) %*% diag(1/sqrt(v))
    if (method=="ZCA-cor") W = G %*% W
  }

  return (W)
}

# whiten data using empirical covariance matrix
whiten = function(X, center=FALSE, method=c("ZCA", "PCA", "Cholesky", "ZCA-cor", "PCA-cor"))
{
    method = match.arg(method)
    S = cov(X)
    W = whiteningMatrix(S, method=method)
    Z = tcrossprod(X, W) # whitened data

    if(center) Z = sweep(Z, 2, colMeans(Z))

    return(Z)
}

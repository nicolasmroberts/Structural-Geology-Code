


# Copyright 2016-2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### GEODESIC REGRESSION ###

# Like oriRescaledGeodesicRegression, but does not rescale its x-data ahead of time.
oriGeodesicRegression <- function(xs, rs, group, numSteps=100) {
  # Let Q be the R whose x is closest to zero.
  q <- rs[[which.min(as.numeric(xs)^2)]]
  # Define the function to be minimized.
  n <- length(rs)
  e <- function(mw) {
    b <- rotExp(rotAntisymmetricFromVector(mw[4:6])) %*% q
    m <- rotAntisymmetricFromVector(mw[1:3])
    f <- function(i) {
      oriDistance(rs[[i]], rotExp(xs[[i]] * m) %*% b, group)^2
    }
    sum(sapply(1:n, f)) / (2 * n)
  }
  # Find the minimum, using the constant geodesic Q as the seed.
  seed <- c(0, 0, 0, 0, 0, 0)
  solution <- optim(seed, e, hessian=TRUE, control=list(maxit=numSteps))
  # Report diagnostic information.
  eigvals <- eigen(solution$hessian, symmetric=TRUE, only.values=TRUE)$values
  m <- rotAntisymmetricFromVector(solution$par[1:3])
  b <- rotExp(rotAntisymmetricFromVector(solution$par[4:6])) %*% q
  rBar <- oriFrechetMean(rs, group)
  rSq <- 1 - solution$value / oriVariance(rs, rBar, group)
  prediction <- function(x) {
    rotExp(x * m) %*% b
  }
  list(m=m, b=b, error=solution$convergence, minEigenvalue=min(eigvals), rSquared=rSq, prediction=prediction)
}

#' Best-fit geodesic curve in SO(3) / G.
#'
#' Returns the best-fit geodesic R(x) = (exp (x m)) b for the given (x, R) pairs. Minimizes the sum of squared distances in the R-direction only, as in ordinary linear regression --- not in the x-direction.
#' @param xs A vector of real numbers.
#' @param rs A list of rotation matrices of same length as xs.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list consisting of $m (anti-symmetric 3x3 real matrix), $b (rotation matrix), $error (real number), $minEigenvalue (real number), $rSquared (real number), $prediction (R function). error is 0 if and only if the minimization succeeds. If error == 1, try increasing numSteps, to 1000 say. minEigenvalue is the least eigenvalue of the Hessian at the solution; worry about the result if it is non-positive 0. rSquared is the R^2 statistic measuring the amount of variance captured, from none (0) to all (1). prediction takes a real number x as input, and returns the predicted rotation matrix R(x) as output.
oriRescaledGeodesicRegression <- function(xs, rs, group, numSteps=100) {
  x0 <- min(xs)
  x1 <- max(xs)
  regr <- oriGeodesicRegression(scales(xs), rs, group, numSteps)
  mNew <- regr$m / (x1 - x0)
  bNew <- rotExp(regr$m * -x0 / (x1 - x0)) %*% regr$b
  prediction <- function(x) {
    regr$prediction((x - x0) / (x1 - x0))
  }
  list(m=mNew, b=bNew, error=regr$error, minEigenvalue=regr$minEigenvalue, rSquared=regr$rSquared, prediction=prediction)
}

#' Permutation test for significance of a geodesic regression.
#' 
#' Obsolete. Use permutedRSquareds instead. Returns up to numPerms R^2 values. May be fewer than numPerms, because of some regressions failing. Let n be the dimension of this vector and g the number of R^2 values greater than the R^2 for the original regression of the data. Let p = g / n. Small values of p indicate that the dependency detected by the regression is meaningful. Uses an internal rescaling, just like oriGeodesicRegression.
#' @param xs A vector of real numbers.
#' @param rs A list of rotation matrices, of same length as xs.
#' @param numPerms A real number (positive integer). The number of permutations to try.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A vector of real numbers, of dimension at most numPerms.
oriGeodesicRegressionPermutations <- function(xs, rs, numPerms, group, numSteps=10000) {
  ys <- scales(xs)
  f <- function(i) {
    print(i / numPerms)
    regr <- oriNativeGeodesicRegression(sample(ys, size=length(ys)), rs, group, numSteps)
    c(regr$error, regr$minEigenvalue, regr$rSquared)
  }
  perms <- sapply(1:numPerms, f)
  perms[3,][perms[1,] == 0 & perms[2,] > 0]
}



### KERNEL REGRESSION ###

#' Kernel regression to fit an orientation R(x) to (x, R) data.
#' 
#' Warning: Sometime in 2017 this function stopped working reliably. I think that the fix is to rescale the xs. Meanwhile, in theory: This function interpolates/extrapolates an orientation R for a given x-value, based on a given set of (x, R) data. The chosen bandwidth h may have a substantial effect on the results. See oriBandwidthForKernelRegression.
#' @param x A real number. The x-value at which to predict the orientation R.
#' @param xs A vector of real numbers. The x-values at which R(x) is known.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param h A real number (positive). The bandwidth. The kernel k is scaled by h, meaning that k is changed to function(x) {k(x / h) / h}.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list consisting of $r (the rotation matrix R(x)), $error (which should be 0; if it is 1, then increase numSteps), and $minEigenvalue (which should be positive; if not, then the minimization has failed).
oriKernelRegression <- function(x, xs, rs, h, group, k=dnorm, numSteps=100) {
  # Let Q be the Ri whose xi is closest to x.
  q <- rs[[which.min((as.numeric(xs) - x)^2)]]
  # Define the function to be minimized.
  kh <- function(x) {k(x / h) / h}
  e <- function(w) {
    r <- rotExp(rotAntisymmetricFromVector(w)) %*% q
    f <- function(i) {
      kh(x - xs[[i]]) * oriDistance(rs[[i]], r, group)^2
    }
    denom <- sum(sapply(xs, function(y) kh(x - y)))
    if (denom == 0)
      1000000
    else {
      numer <- sum(sapply(1:length(xs), f))
      numer / denom
    }
  }
  # Minimize the function with respect to R, using Q as the seed.
  seed <- c(0, 0, 0)
  solution <- optim(seed, e, hessian=TRUE, control=list(maxit=numSteps))
  # Report diagnostic information.
  eigvals <- eigen(solution$hessian, symmetric=TRUE, only.values=TRUE)$values
  r <- rotExp(rotAntisymmetricFromVector(solution$par)) %*% q
  list(r=r, error=solution$convergence, minEigenvalue=min(eigvals))
}

#' Cross-validation algorithm for choosing the bandwidth for kernel regression.
#' 
#' Warning: Sometime in 2017 this function stopped working reliably. See oriKernelRegression.
#' @param xs A vector of real numbers. The x-values at which R(x) is known. Assumed to have length >= 3.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @param lowerBound A real number (positive). A lower bound on the bandwidth.
#' @param upperBound A real number (positive). An upper bound on the bandwidth.
#' @return A real number (positive). The bandwidth to use in kernel regression.
oriBandwidthForKernelRegression <- function(xs, rs, group, k=dnorm, numSteps=100, lowerBound=0.00001, upperBound=pi) {
  # Build jackknifed lists ahead of time.
  n <- length(xs)
  xjs <- listOmitting(xs)
  rjs <- listOmitting(rs)
  # Define the function to minimize.
  g <- function(h) {
    rjhs <- lapply(1:n, function(j) oriKernelRegression(xs[[j]], xjs[[j]], rjs[[j]], h=h, group=group, k=k, numSteps=numSteps))
    sum(sapply(1:n, function(j) oriDistance(rs[[j]], rjhs[[j]]$r, group)^2))
  }
  # Minimize the function with respect to log h.
  solution <- optimize(g, lower=lowerBound, upper=upperBound)
  solution$minimum
}

#' R^2 statistic for kernel regression.
#' 
#' Warning: Sometime in 2017 this function stopped working reliably. See oriKernelRegression.
#' @param xs A vector of real numbers. The x-values at which R(x) is known.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param h A real number (positive). The bandwidth. The kernel k is scaled by h, meaning that k is changed to function(x) {k(x / h) / h}.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A real number (usually between 0 and 1). The R^2 statistic.
oriRsquaredForKernelRegression <- function(xs, rs, h, group, k=dnorm, numSteps=100) {
  n <- length(xs)
  qs <- lapply(xs, function(x) oriKernelRegression(x, xs, rs, h=h, group=group, k=k, numSteps=numSteps)$r)
  e <- sum(sapply(1:n, function(i) oriDistance(rs[[i]], qs[[i]], group)^2)) / (2 * n)
  1 - e / oriVariance(rs, oriFrechetMean(rs, group), group)
}



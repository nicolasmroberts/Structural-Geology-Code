


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### GEODESIC REGRESSION ###

#' Regression to fit a geodesic curve to data (x, R), with no pre-processing of the independent variable.
#'
#' Returns the best-fit geodesic R(x) = (exp (x m)) b for the given (x, r) pairs. Minimizes the sum of squared distances in the r-direction only, as in ordinary linear regression. Usually you want to try rotRescaledGeodesicRegression first.
#' @param xs A list of real numbers.
#' @param rs A list of rotation matrices of same length as xs.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list consisting of $m (anti-symmetric 3x3 real matrix), $b (rotation matrix), $error (real number), $minEigenvalue (real number), $rSquared (real number), $prediction (R function). $error is 0 if and only if the minimization succeeds. If $error == 1, try increasing steps, to 1000 say. $minEigenvalue is the least eigenvalue of the Hessian at the solution; worry about the result if it is non-positive. $rSquared is the R^2 statistic measuring the amount of variance captured, from none (0) to all (1). $prediction takes an x as input and returns the predicted rotation R(x).
rotGeodesicRegression <- function(xs, rs, numSteps=100) {
  # Let Q be the R whose x is closest to zero.
  q <- rs[[which.min(as.numeric(xs)^2)]]
  # Define the function to be minimized.
  n <- length(rs)
  e <- function(mw) {
    b <- rotExp(rotAntisymmetricFromVector(mw[4:6])) %*% q
    m <- rotAntisymmetricFromVector(mw[1:3])
    f <- function(i) {
      rotDistance(rs[[i]], rotExp(xs[[i]] * m) %*% b)^2
    }
    sum(sapply(1:n, f)) / (2 * n)
  }
  # Find the minimum, using the constant geodesic Q as the seed.
  seed <- c(0, 0, 0, 0, 0, 0)
  solution <- optim(seed, e, lower=c(-pi, -pi, -pi, 0), upper=c(pi, pi, pi, pi), hessian=TRUE, control=list(maxit=numSteps))
  # Report diagnostic information.
  eigvals <- eigen(solution$hessian, symmetric=TRUE, only.values=TRUE)$values
  m <- rotAntisymmetricFromVector(solution$par[1:3])
  b <- rotExp(rotAntisymmetricFromVector(solution$par[4:6])) %*% q
  rSq <- 1 - solution$value / rotVariance(rs, rotFrechetMean(rs))
  prediction <- function(x) {
    rotExp(x * m) %*% b
  }
  list(b=b, m=m, error=solution$convergence, minEigenvalue=min(eigvals), rSquared=rSq, prediction=prediction)
}

#' Geodesic curve regression, with a convenient behind-the-scenes rescaling.
#' 
#' In theory, this function is equivalent to rotGeodesicRegression. In practice, it uses a behind-the-scenes rescaling of the xs onto the interval [0, 1], which I find improves the performance of the optimization. You should use this version of the function, unless you have some specific reason to use rotGeodesicRegression.
#' @param xs A list of real numbers.
#' @param rs A list of rotation matrices of same length as xs.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list consisting of $m (anti-symmetric 3x3 real matrix), $b (rotation matrix), $error (real number), $minEigenvalue (real number), $rSquared (real number), $prediction (R function). $error is 0 if and only if the minimization succeeds. If $error == 1, try increasing steps, to 1000 say. $minEigenvalue is the least eigenvalue of the Hessian at the solution; worry about the result if it is non-positive. $rSquared is the R^2 statistic measuring the amount of variance captured, from none (0) to all (1). $prediction takes an x as input and returns the predicted rotation R(x).
rotRescaledGeodesicRegression <- function(xs, rs, numSteps=100) {
  x0 <- min(xs)
  x1 <- max(xs)
  regr <- rotGeodesicRegression(scales(xs), rs, numSteps)
  mNew <- regr$m / (x1 - x0)
  bNew <- rotExp(regr$m * -x0 / (x1 - x0)) %*% regr$b
  prediction <- function(x) {
    regr$prediction((x - x0) / (x1 - x0))
  }
  list(m=mNew, b=bNew, error=regr$error, minEigenvalue=regr$minEigenvalue, rSquared=regr$rSquared, prediction=prediction)
}

#' Permutation test for significance of a geodesic regression.
#' 
#' Obsolete. Use permutedRSquareds instead. Returns up to numPerms R^2 values. May be fewer than numPerms, because of some regressions failing. Let n be the dimension of this vector and g the number of R^2 values greater than the R^2 for the original regression of the data. Let p = g / n. Small values of p indicate that the dependency detected by the regression is meaningful. Uses an internal rescaling, just like rotGeodesicRegression.
#' @param xs A vector of real numbers.
#' @param rs A list of rotation matrices, of same length as xs.
#' @param numPerms A real number (positive integer). The number of permutations to try.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A vector of real numbers, of dimension at most numPerms.
rotGeodesicRegressionPermutations <- function(xs, rs, numPerms, numSteps=10000) {
  ys <- scales(xs)
  f <- function(i) {
    print(i / numPerms)
    regr <- rotNativeGeodesicRegression(sample(ys, size=length(ys)), rs, numSteps)
    c(regr$error, regr$minEigenvalue, regr$rSquared)
  }
  perms <- sapply(1:numPerms, f)
  perms[3,][perms[1,] == 0 & perms[2,] > 0]
}



### KERNEL REGRESSION ###

#' Kernel regression to fit a rotation R(x) to (x, R) data, with no pre-processing of the independent variable.
#' 
#' Warning: Sometime in 2017 this function stopped working reliably. I think that the solution is to rescale the xs. Stay tuned. Theoretically identical to rotRescaledKernelRegression, but doesn't rescale the xs to the interval [0, 1].
#' @param x A real number. The x-value at which to predict the rotation R.
#' @param xs A vector of real numbers. The x-values at which R(x) is known.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param h A real number (positive). The bandwidth. The kernel k is scaled by h, meaning that k is changed to function(x) {k(x / h) / h}.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list consisting of $r (the rotation matrix R(x)), $error (which should be 0; if it is 1, then increase numSteps), and $minEigenvalue (which should be positive; if not, then the minimization has failed).
rotKernelRegression <- function(x, xs, rs, h, k=dnorm, numSteps=100) {
  # Let Q be the Ri whose xi is closest to x.
  q <- rs[[which.min((xs - x)^2)]]
  # Define the function to be minimized.
  kh <- function(x) {k(x / h) / h}
  e <- function(w) {
    r <- rotExp(rotAntisymmetricFromVector(w)) %*% q
    f <- function(i) {
      kh(x - xs[[i]]) * rotDistance(rs[[i]], r)^2
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
#' Warning: Sometime in 2017 this function stopped working reliably. See rotKernelRegression. See Davis et al. (2010).
#' @param xs A vector of real numbers. The x-values at which R(x) is known. Assumed to have length >= 3.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A real number (positive). The bandwidth to use in kernel regression.
rotBandwidthForKernelRegression <- function(xs, rs, k=dnorm, numSteps=100) {
  # Build jackknifed lists ahead of time.
  n <- length(xs)
  xjs <- listOmitting(xs)
  rjs <- listOmitting(rs)
  # Define the function to minimize.
  g <- function(h) {
    rjhs <- lapply(1:n, function(j) rotKernelRegression(xs[[j]], xjs[[j]], rjs[[j]], h, k, numSteps))
    sum(sapply(1:n, function(j) rotDistance(rs[[j]], rjhs[[j]]$r)^2))
  }
  # Minimize the function with respect to h.
  solution <- optimize(g, interval=c(0, pi))
  solution$minimum
}

# What is the point of this? It might be orphaned.
rotKernelRegressionBandwidthMisfit <- function(xs, rs, h, k=dnorm, numSteps=100) {
  n <- length(xs)
  xjs <- listOmitting(xs)
  rjs <- listOmitting(rs)
  rjhs <- lapply(1:n, function(j) rotKernelRegression(xs[[j]], xjs[[j]], rjs[[j]], h, k, numSteps))
  sum(sapply(1:n, function(j) rotDistance(rs[[j]], rjhs[[j]]$r)^2))
}

#' R^2 statistic for kernel regression.
#' 
#' Warning: Sometime in 2017 this function stopped working reliably. See rotKernelRegression. As always, make sure the errors are 0 and the minEigenvalues are positive.
#' @param xs A vector of real numbers. The x-values at which R(x) is known.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param h A real number (positive). The bandwidth. The kernel k is scaled by h, meaning that k is changed to function(x) {k(x / h) / h}.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list with elements $rSquared, $errors, $minEigenvalues.
rotRsquaredForKernelRegression <- function(xs, rs, h, k=dnorm, numSteps=100) {
  n <- length(xs)
  regrs <- lapply(xs, function(x) rotKernelRegression(x, xs, rs, h, k, numSteps))
  e <- sum(sapply(1:n, function(i) rotDistance(rs[[i]], regrs[[i]]$r)^2)) / (2 * n)
  rSquared <- 1 - e / rotVariance(rs, rotFrechetMean(rs))
  errors <- sapply(regrs, function(regr) regr$error)
  minEigenvalues <- sapply(regrs, function(regr) regr$minEigenvalue)
  list(rSquared=rSquared, errors=errors, minEigenvalues=minEigenvalues)
}

#' Permutation test for significance of a kernel regression.
#' 
#' Warning: Sometime in 2017 this function stopped working reliably. See rotKernelRegression. Returns up to numPerms R^2 values. May be fewer than numPerms, because of some regressions failing. Let n be the dimension of this vector and g the number of R^2 values greater than the R^2 for the original regression of the data. Let p = g / n. Small values of p indicate that the dependency detected by the regression is meaningful.
#' @param xs A vector of real numbers. The x-values at which R(x) is known.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param h A real number (positive). The bandwidth. The kernel k is scaled by h, meaning that k is changed to function(x) {k(x / h) / h}.
#' @param numPerms A real number (positive integer). The number of permutations to try.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A vector of length up to numPerms, containing R^2 values for successful permutations.
rotKernelRegressionPermutations <- function(xs, rs, h, numPerms, k=dnorm, numSteps=10000) {
  ys <- scales(xs)
  f <- function(i) {
    print(i / numPerms)
    rSquaredErrorsMins <- rotRsquaredForKernelRegression(sample(ys, size=length(ys)), rs, h, numSteps=numSteps)
    error <- max(abs(rSquaredErrorsMins$errors))
    minEigenvalue <- min(rSquaredErrorsMins$minEigenvalues)
    c(error, minEigenvalue, rSquaredErrorsMins$rSquared)
  }
  perms <- sapply(1:numPerms, f)
  perms[3,][perms[1,] == 0 & perms[2,] > 0]
}

#' Kernel regression to fit a rotation R(x) to (x, R) data, with no pre-processing of the independent variable.
#' 
#' Warning: Sometime in 2017 this function stopped working reliably. See rotKernelRegression. This function interpolates/extrapolates a rotation R for a given x-value, based on a given set of (x, R) data. The chosen bandwidth h may have a substantial effect on the results. See rotBandwidthForKernelRegression. If this function doesn't work as you expect, then try rotKernelRegression.
#' @param x A real number. The x-value at which to predict the rotation R.
#' @param xs A vector of real numbers. The x-values at which R(x) is known.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param h A real number (positive). The bandwidth. The kernel k is scaled by h, meaning that k is changed to function(x) {k(x / h) / h}.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list consisting of $r (the rotation matrix R(x)), $error (which should be 0; if it is 1, then increase numSteps), and $minEigenvalue (which should be positive; if not, then the minimization has failed).
rotRescaledKernelRegression <- function(x, xs, rs, h, k=dnorm, numSteps=100) {
  x0 <- min(xs)
  x1 <- max(xs)
  rotKernelRegression((x - x0) / (x1 - x0), scales(xs), rs, h, k, numSteps)
}



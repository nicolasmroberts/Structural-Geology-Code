


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### MATRIX FISHER DISTRIBUTION ###

#' Maximum likelihood estimation of the matrix Fisher parameters M and K.
#'
#' The method proceeds by minimizing a certain function using numerical methods. Things can go wrong. The results include two kinds of diagnostic information. Warning: Don't be surprised if this doesn't work very well yet. The seeds need better choosing. Also, we should constrain the optimization to prevent the integrals from getting insanely big.
#' @param rs A list of rotation matrices.
#' @param seeds A list of 3-dimensional real vectors, or NULL. Seeds for the minimization. The details are complicated. If NULL is given, then this function picks 6 of them automatically.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list consisting of $mHat (rotation matrix), $kHat (symmetric 3x3 real matrix), $error (real number), $minEigenvalue (real number). $error is 0 if and only if the minimization succeeds. If $error is 1, then try increasing numSteps. $minEigenvalue is the least eigenvalue of the Hessian at the solution; worry about the result if it is non-positive.
rotFisherMLE <- function(rs, seeds=NULL, numSteps=1000) {
  # Shockingly, the data are boiled down to three numbers for all of the hard stuff.
  udv <- rotSingularValueDecomposition(arithmeticMean(rs))
  d <- diag(udv$d)
  # Here is the function to be minimized.
  f <- function(logP) {
    p <- exp(logP)
    s <- sort(p, decreasing=TRUE)
    f <- function(u) {
      besselI((s[1] - s[2]) * u, 0) * besselI((s[1] + s[2]) * (1 - u), 0) * exp(s[3] * (1 - 2 * u))
    }
    log(integrate(f, 0, 1)$value) - p %*% d
  }
  # Choose seeds based on the first octant of the 24 chambers of Sei et al. (2013, p. 448).
  if (is.null(seeds))
    seeds <- list(c(-1, 0, 1), c(-1, 1, 0), c(0, -1, 1), c(0, 1, -1), c(1, -1, 0), c(1, 0, -1))
  # Try all of the seeds and report the best answer found.
  solutions <- lapply(seeds, function(seed) optim(seed, f, hessian=TRUE, control=list(maxit=numSteps)))
  values <- sapply(solutions, function(solution) solution$value)
  solution <- solutions[[which.min(values)]]
  mHat <- udv$u %*% t(udv$v)
  kHat <- udv$v %*% diag(exp(solution$par)) %*% t(udv$v)
  eigvals <- eigen(solution$hessian, symmetric=TRUE, only.values=TRUE)$values
  list(mHat=mHat, kHat=kHat, error=solution$convergence, minEigenvalue=min(eigvals))
}

#' P-value function for whether the data come from a matrix Fisher distribution with the hypothesized mean.
#'
#' Based on Downs (1972; Eq. 5.7, p = 2 case). Assumes large sample size and tightly clustered sample. Uses the Stiefel manifold version of the matrix Fisher distribution, rather than the SO(3) version (Sei et al. (2013)). Only the left 3x2 submatrices of the given matrices are used. Indeed, the matrices may safely be given as 3x2.
#' @param rs A list of rotation matrices.
#' @return An R function from {rotation matrices} to {real numbers}. For any given hypothesized mean, this function produces the p-value for that hypothesis.
rotDownsInference <- function(rs) {
  n <- length(rs)
  r32s <- lapply(rs, function(a) {a[1:3, 1:2]})
  # hHat = (rBar^T rBar)^(1 / 2).
  rBar <- arithmeticMean(r32s)
  eig <- eigen(crossprod(rBar, rBar), symmetric=TRUE)
  hHat <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  top <- 2.0 * (diag(c(1.0, 1.0)) - hHat)
  top <- sqrt(tr(crossprod(top, top)))
  f <- function(r) {
    r32 <- r[1:3, 1:2]
    bottom <- diag(c(2, 2)) - crossprod(rBar, r32) - crossprod(r32, rBar)
    xi <- top / sqrt(tr(crossprod(bottom, bottom)))
    q <- (1 / sqrt(xi) - 1) * (n - 5 / 3)
    1 - pf(q, 3, 3 * n - 5)
  }
  f
}

#' One rotation drawn from the matrix Fisher distribution on SO(3).
#'
#' This function is a helper function for rotFisher. Probably you do not want to call this function yourself.
#' @param lambda A 4D real vector. lambda from Kent et al. (2013).
#' @param omega A 4x4 real matrix. Omega from Kent et al. (2013).
#' @param sigma A 4x4 real matrix. Inverse of omega. Sigma from Kent et al. (2013).
#' @param bound A real number. M* from Kent et al. (2013).
#' @return A rotation matrix.
rotFisherHelper <- function(lambda, omega, sigma, bound) {
  x <- rayNormalized(mvrnorm(1, c(0, 0, 0, 0), sigma))
  w <- runif(1, 0, 1)
  # Repeat until w < f*(x) / (M* g*(x)).
  while (bound * w >= exp(-x %*% lambda %*% x) * (x %*% omega %*% x)^2) {
    x <- rayNormalized(mvrnorm(1, c(0, 0, 0, 0), sigma))
    w <- runif(1, 0, 1)
  }
  rotMatrixFromQuaternion(x)
}

#' Sampling from the matrix Fisher distribution on SO(3).
#'
#' Follows Kent et al. (2013), 'A new method to simulate the Bingham and related distributions in directional data analysis with applications'.
#' @param m A rotation matrix. The center of the distribution.
#' @param k A 3x3 real matrix (symmetric, positive-definite). The concentration of the distribution.
#' @param n A real number (positive integer) or NULL.
#' @return Either a rotation matrix (if n is NULL) or a list of n rotation matrices (if n is a number).
rotFisher <- function(m, k, n=NULL) {
  if (is.null(n))
    rotFisher(m, k, 1)[[1]]
  else {
    # Compute the Lambdas for the Bingham distribution.
    eig <- eigen(k, symmetric=TRUE)
    d <- eig$values
    l1 <- 2 * (d[[2]] + d[[3]])
    l2 <- 2 * (d[[1]] + d[[3]])
    l3 <- 2 * (d[[1]] + d[[2]])
    lams <- c(0, l1, l2, l3)
    # Compute b as the greatest real root of a certain quartic polynomial.
    a0 <- 8 * l1 * l2 * l3
    a1 <- 8 * (l1 * l2 + l1 * l3 + l2 * l3 - l1 * l2 * l3)
    a2 <- 6 * (l1 + l2 + l3) - 4 * (l1 * l2 + l1 * l3 + l2 * l3)
    a3 <- 4 - 2 * (l1 + l2 + l3)
    a4 <- -1
    roots <- polyroot(c(a0, a1, a2, a3, a4))
    roots <- Filter(function (x) (Im(x)^2 < epsilon), roots)
    roots <- sapply(roots, function(x) Re(x))
    b <- max(roots)
    # Compute bound M*, Omega for ACG, and Sigma = Omega^-1 for normal.
    bound <- exp((b - 4) / 2) * (4 / b)^2
    oms <- c(1, 1, 1, 1) + 2 * lams / b
    sigs <- 1 / oms
    # Generate the random sample.
    v <- eig$vectors
    rs <- replicate(n, rotFisherHelper(diag(lams), diag(oms), diag(sigs), bound), simplify=FALSE)
    lapply(rs, function(r) (m %*% v %*% r %*% t(v)))
  }
}

#' Inference about the population mean, based on parametric bootstrapping with the Fisher distribution.
#' 
#' Identical to rotBootstrapInference, but draws its bootstrap samples from the MLE matrix Fisher distribution for the data, rather than from the data set itself.
#' @param rs A list of rotation matrices.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @param func An R function. Either rotFrechetMean or rotProjectedMean (or something equivalent).
#' @return A list consisting of $bootstraps (a list of rotation matrices) and everything returned by rotMahalanobisInference.
rotFisherBootstrapInference <- function(rs, numBoots, func=rotFrechetMean) {
  mle <- rotFisherMLE(rs, numSteps=10000)
  boots <- replicate(numBoots, func(rotFisher(mle$mHat, mle$kHat, length(rs)), ...), simplify=FALSE)
  bootMean <- func(boots, ...)
  infer <- rotMahalanobisInference(boots, bootMean)
  infer$bootstraps <- boots
  infer
}

#' Two-sample inference about the difference in population means, based on parametric bootstrapping using the Fisher distribution.
#' 
#' Identical to rotTwoSampleBootstrapInference, but draws its bootstrap samples from matrix Fisher distributions MLE-fit to the data, rather than from the data set itself.
#' @param firsts A list of rotation matrices.
#' @param seconds A list of rotation matrices.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @param func An R function. Either rotFrechetMean or rotProjectedMean (or something equivalent).
#' @return A list consisting of $bootstraps (a list of rotation matrices) and everything returned by rotMahalanobisInference.
rotTwoSampleFisherBootstrapInference <- function(firsts, seconds, numBoots, func=rotFrechetMean) {
  firstMLE <- rotFisherMLE(firsts, numSteps=10000)
  secondMLE <- rotFisherMLE(seconds, numSteps=10000)
  f <- function() {
    firstMean <- func(rotationFisher(firstMLE$mHat, firstMLE$kHat, length(firsts)))
    secondMean <- func(rotationFisher(secondMLE$mHat, secondMLE$kHat, length(seconds)))
    secondMean %*% t(firstMean)
  }
  boots <- replicate(numBoots, f(), simplify=FALSE)
  bootMean <- func(boots, ...)
  infer <- rotMahalanobisInference(boots, bootMean)
  infer$bootstraps <- boots
  infer
}



### ISOTROPIC MATRIX FISHER DISTRIBUTION ###

#' Sampling from the isotropic matrix Fisher distribution.
#'
#' @param s A rotation matrix. The center of the distribution.
#' @param kappa A real number (positive). The concentration of the distribution.
#' @param n A real number (positive integer) or NULL.
#' @return Either a rotation matrix (if n is NULL) or a list of n rotation matrices (if n is a number).
rotIsotropicFisher <- function(s, kappa, n=NULL) {
  rotFisher(s, diag((kappa^2 / 2) * c(1, 1, 1)), n)
}

#' Jeffreys prior for the concentration parameter of the isotropic matrix Fisher distribution on SO(3).
#'
#' See Bingham et al. (2010, Section 3.2). To avoid a notational conflict, rename their kappa to lambda. The relationship between lambda and the kappa of Qiu et al. (2013) is lambda == kappa^2 / 2. Then also let eta == -log kappa, so that kappa = exp(-eta). We work in terms of eta.
rotIsotropicFisherJeffreysPrior <- function(negLogKappa) {
  if (177 < negLogKappa)
    0
  else if (negLogKappa <= -3.433049)
    sqrt(6)
  else if (negLogKappa < -2.94) {
    # Linearly interpolate the function between -3.433049 and -2.94.
    y1 <- 2.451246 # == rotIsotropicFisherJeffreysPrior(-2.93)
    y0 <- 2.451212 # == rotIsotropicFisherJeffreysPrior(-2.94)
    m <- (y1 - y0) /  (-2.93 - -2.94)
    m * (negLogKappa + 2.94) + y0
  } else {
    #lam <- exp(-2 * negLogKappa) / 2
    #i0 <- besselI(2 * lam, 0)
    #i1 <- besselI(2 * lam, 1)
    #numer <- i0^2 * 2 / lam - i0 * i1 * 2 / lam^2 + i1^2 * (1 / lam^2 - 2 / lam)
    #denom <- (i0 - i1)^2
    #sqrt(numer / denom) * 2 * lam
    kSq <- exp(-2 * negLogKappa)
    i0 <- besselI(kSq, 0)
    i1 <- besselI(kSq, 1)
    numer <- (kSq - 1) * i0^2 - kSq * i1^2
    denom <- (i0 - i1)^2
    jeff <- 2 * sqrt(1 + numer / denom)
    jeff
  }
}

# This code reproduces the Fisher panel of Fig. 2 of Qiu et al. (2013).
#etas <- seq(from=-2, to=5, by=0.01)
#priors <- sapply(etas, rotIsotropicFisherJeffreysPrior)
#plot(x=etas, y=(priors / sqrt(6)))

# This code checks the behavior around 177. Looks good.
#etas <- seq(from=170, to=180, by=0.01)
#priors <- sapply(etas, rotIsotropicFisherJeffreysPrior)
#plot(x=etas, y=priors)

# This code checks the behavior around -3.433 and -2.94. Not smooth, but okay.
#etas <- seq(from=-4, to=-2, by=0.01)
#priors <- sapply(etas, rotIsotropicFisherJeffreysPrior)
#plot(x=etas, y=priors)



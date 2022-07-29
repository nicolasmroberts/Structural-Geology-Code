


# Copyright 2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### INFRASTRUCTURE ###

# Returns C such that C exp(-u^T A u) integrates to 1 over the hemisphere? Accuracy is doubtful.
lineBinghamSaddlepointConstant <- function(a) {
  vals <- eigen(a, symmetric=TRUE, only.values=TRUE)$values
  const <- fb.saddle(c(0, 0, 0), vals)
  1 / const[[3]]
}

# Returns C such that C exp(-u^T A u) integrates to 1 over the hemisphere? Accuracy is unclear.
lineBinghamIntegratedConstant <- function(a, numNonAdapt=5) {
  rayInt <- rayIntegral(function(u) exp(-u %*% a %*% u), numNonAdapt)
  1 / (0.5 * rayInt)
}

# Test. The two methods seem to disagree in a systematic way.
lineBinghamConstantTest <- function(sd) {
  a <- rnorm(5, mean=0, sd=sd)
  a <- rbind(c(a[[1]], a[[2]], a[[3]]), c(a[[2]], a[[4]], a[[5]]), c(a[[3]], a[[5]], -a[[1]] - a[[4]]))
  c(1 / lineBinghamSaddlepointConstant(a), 1 / lineBinghamIntegratedConstant(a))
}
#plot(t(replicate(100, lineBinghamConstantTest(2))))



### MISCELLANY ###

#' Simulation from the Bingham distribution.
#' 
#' Uses the function rbingham in package Directional. In that convention, the Bingham probability density is proportional to exp(-x^T A x), not exp(x^T A x). The mean is the eigenvector of A with least eigenvalue. The main direction of the dispersion is toward the eigenvector of A with the intermediate eigenvalue. The pole to the dispersion is the eigenvector with greatest eigenvalue.
#' @param a A symmetric real 3x3 matrix.
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single line. If n is a positive integer, then a list of n lines.
lineBingham <- function(a, n=NULL) {
  if (is.null(n))
    lower(as.numeric(rbingham(1, a)))
  else {
    us <- rbingham(n, a)
    if (n == 1)
      list(lower(as.numeric(us)))
    else
      lapply(1:nrow(us), function(i) lower(us[i,]))
  }
}

#' Maximum likelihood estimation of the Bingham distribution parameters.
#' 
#' Compared to lineTauxeBinghamMLE, this is slow but apparently more accurate. Uses a numerical integration (to compute the normalization constant) inside an numerical optimization (to maximize the likelihood). The Bingham probability density is proportional to exp(-x^T A x), not exp(x^T A x). I prefer to normalize A so that tr A == 0.
#' @param us A list of lines.
#' @param weights A vector of real numbers (non-negative), of length equal to us. They need not sum to 1; the function automatically normalizes them to do so.
#' @param numNonAdapt A real number (non-negative integer). The number of refinements to use in the numerical integration. Each increment of numNonAdapt increases time and memory requirements by a factor of 4.
#' @param numSteps. A real number (positive integer). The number of steps to use in the numerical optimization. If the output $error is non-zero, then try increasing numSteps.
#' @return A list with members $a (symmetric 3x3 real matrix A), $values (a real 3D vector; the eigenvalues of A; sum to zero), $vectors (a rotation matrix; the eigenvectors of A are the columns), $error (integer; increase numSteps if $error != 0), and $minEigenvalue (the minimum eigenvalue of the Hessian at the putative optimum; worry if this is not positive).
lineBinghamMLE <- function(us, weights=replicate(length(us), 1), numNonAdapt=5, numSteps=1000) {
  # Scale the weights so that they sum to 1.
  n <- length(us)
  ws <- weights / sum(weights)
  # In eigen, the eigenvectors are descending and the eigenvectors are unit.
  ts <- lapply(1:n, function(i) {ws[[i]] * outer(us[[i]], us[[i]])})
  tBar <- Reduce("+", ts)
  eig <- eigen(tBar, symmetric=TRUE)
  # Let's assume that the MLE A will have the same eigenvectors as T-bar, as in the unweighted case.
  q <- eig$vectors
  vs <- lapply(us, function(u) as.numeric(t(q) %*% u))
  # The MLE will have eigenvalues a1, a2, -a1 - a2. Define f to be the negative log-likelihood.
  f <- function(a1a2) {
    diagonal <- diag(c(a1a2[[1]], a1a2[[2]], -a1a2[[1]] - a1a2[[2]]))
    const <- lineBinghamIntegratedConstant(diagonal, numNonAdapt)
    #-log(const) + sum(sapply(vs, function(v) {v %*% diagonal %*% v}))
    -log(const) + sum(sapply(1:n, function(i) {ws[[i]] * vs[[i]] %*% diagonal %*% vs[[i]]}))
  }
  # Numerically minimize f to obtain the maximum likelihood estimate.
  seed <- c(0, 0)
  sol <- optim(seed, f, hessian=TRUE, control=list(maxit=numSteps))
  values <- c(sol$par[[1]], sol$par[[2]], -sol$par[[1]] - sol$par[[2]])
  a <- q %*% diag(values) %*% t(q)
  # Report diagnostic information.
  eigvals <- eigen(sol$hessian, symmetric=TRUE, only.values=TRUE)$values
  list(error=sol$convergence, minEigenvalue=min(eigvals), a=a, values=values, vectors=q)
}

# Test. lineBinghamMLE works well when I use my integrated constant, and terribly when I use the saddlepoint-approximation constant from the package Directional.
lineBinghamMLETest <- function(n, sd) {
  a <- rnorm(5, mean=0, sd=sd)
  a <- rbind(c(a[[1]], a[[2]], a[[3]]), c(a[[2]], a[[4]], a[[5]]), c(a[[3]], a[[5]], -a[[1]] - a[[4]]))
  us <- lineBingham(a, n)
  mleDavis <- lineBinghamMLE(us, numSteps=10000)
  print(c(mleDavis$error, mleDavis$minEigenvalue))
  mleTauxe <- lineTauxeBinghamMLE(us)
  print(eigen(a, symmetric=TRUE, only.values=TRUE)$values)
  print(eigen(mleDavis$a, symmetric=TRUE, only.values=TRUE)$values)
  print(mleTauxe$values)
  print(a)
  print(mleDavis$a)
  print(mleTauxe$a)
}
#lineBinghamMLETest(1000, 4)



### TAUXE (2010) ###

# Bingham MLE lookup table from Appendix C of Tauxe's books.
lineBinghamK1Table <- cbind(
  c(-25.55, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.56, -13.11, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.58, -13.14, -9.043, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.6, -13.16, -9.065, -7.035, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.62, -13.18, -9.080, -7.042, -5.797, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.63, -13.19, -9.087, -7.041, -5.789, -4.917, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.64, -13.20, -9.087, -7.033, -5.773, -4.896, -4.231, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.65, -13.20, -9.081, -7.019, -5.752, -4.868, -4.198, -3.659, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.65, -13.19, -9.068, -6.999, -5.726, -4.836, -4.160, -3.616, -3.160, NA, NA, NA, NA, NA, NA, NA),
  c(-25.64, -13.18, -9.05, -6.974, -5.694, -4.799, -4.118, -3.570, -3.109, -2.709, NA, NA, NA, NA, NA, NA),
  c(-25.63, -13.17, -9.027, -6.944, -5.658, -4.757, -4.071, -3.518, -3.053, -2.649, -2.289, NA, NA, NA, NA, NA),
  c(-25.61, -23.14, -8.999, -6.910, -5.618, -4.711, -4.021, -3.463, -2.993, -2.584, -2.220, -1.888, NA, NA, NA, NA),
  c(-25.59, -13.12, -8.966, -6.870, -5.573, -4.661, -3.965, -3.403, -2.928, -2.515, -2.146, -1.809, -1.497, NA, NA, NA),
  c(-25.57, -13.09, -8.928, -6.827, -5.523, -4.606, -3.906, -3.338, -2.859, -2.441, -2.066, -1.724, -1.406, -1.106, NA, NA),
  c(-25.54, -13.05, -8.886, -6.778, -5.469, -4.547, -3.842, -3.269, -2.785, -2.361, -1.981, -1.634, -1.309, -1.002, -0.708, NA),
  c(-25.50, -13.01, -8.839, -6.725, -5.411, -4.484, -3.773, -3.195, -2.706, -2.277, -1.891, -1.537, -1.206, -0.891, -0.588, -0.292),
  c(-25.46, -12.96, -8.788, -6.668, -5.348, -4.415, -3.699, -3.116, -2.621, -2.186, -1.794, -1.433, -1.094, -0.771, -0.459, -0.152),
  c(-25.42, -12.91, -8.731, -6.606, -5.280, -4.342, -3.620, -3.032, -2.531, -2.089, -1.690, -1.322, -0.974, -0.642, NA, NA),
  c(-25.37, -12.86, -8.670, -6.539, -5.207, -4.263, -3.536, -2.941, -2.434, -1.986, -1.579, -1.202, NA, NA, NA, NA),
  c(-25.31, -12.80, -8.604, -6.466, -5.126, -4.179, -3.446, -2.845, -2.330, -1.874, NA, NA, NA, NA, NA, NA),
  c(-25.5, -12.73, -8.532, -6.388, -5.045, -4.089, -3.349, -2.741, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.19, -12.66, -8.454, -6.305, -4.955, -3.992, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.12, -12.58, -8.371, -6.215, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))

# Bingham MLE lookup table from Appendix C of Tauxe's books.
lineBinghamK2Table <- cbind(
  c(-25.55, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-13.09, -13.11, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-8.996, -9.019, -9.043, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-6.977, -6.999, -7.020, -7.035, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-5.760, -5.777, -5.791, -5.798, -5.797, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-4.923, -4.934, -4.941, -4.941, -4.933, -4.917, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-4.295, -4.301, -4.301, -4.294, -4.279, -4.258, -4.231, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-3.796, -3.796, -3.790, -3.777, -3.756, -3.729, -3.697, -3.659, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-3.381, -3.375, -3.363, -3.345, -3.319, -3.287, -3.249, -3.207, -3.160, NA, NA, NA, NA, NA, NA, NA),
  c(-3.025, -3.014, -2.997, -2.973, -2.942, -2.905, -2.863, -2.816, -2.765, -2.709, NA, NA, NA, NA, NA, NA),
  c(-2.712, -2.695, -2.673, -2.644, -2.609, -2.568, -2.521, -2.470, -2.414, -2.354, -2.289, NA, NA, NA, NA, NA),
  c(-2.431, -2.410, -2.382, -2.349, -2.309, -2.263, -2.212, -2.157, -2.097, -2.032, -1.963, -1.888, NA, NA, NA, NA),
  c(-2.175, -2.149, -2.117, -2.078, -2.034, -1.984, -1.929, -1.869, -1.805, -1.735, -1.661, -1.582, -1.497, NA, NA, NA),
  c(-1.939, -1.908, -1.871, -1.828, -1.779, -1.725, -1.665, -1.601, -1.532, -1.458, -1.378, -1.294, -1.203, -1.106, NA, NA),
  c(-1.718, -1.682, -1.641, -1.596, -1.540, -1.481, -1.417, -1.348, -1.274, -1.195, -1.110, -1.020, -0.923, -0.819, -0.708, NA),
  c(-1.510, -1.470, -1.423, -1.371, -1.313, -1.250, -1.181, -1.108, -1.028, -0.944, -0.853, -0.756, -0.653, -0.541, -0.421, -0.292), 
  c(-1.312, -1.267, -1.216, -1.159, -1.096, -1.028, -0.955, -0.876, -0.791, -0.701, -0.604, -0.500, -0.389, -0.269, -0.140, 0.000),
  c(-1.123, -1.073, -1.017, -9.555, -0.887, -0.814, -0.736, -0.651, -0.561, -0.464, -0.360, -0.249, -0.129, 0.000, NA, NA),
  c(-0.940, -0.885, -0.824, -0.757, -0.684, -0.606, -0.522, -0.432, -0.335, -0.231, -0.120, 0.000, NA, NA, NA, NA),
  c(-0.762, -0.702, -0.636, -0.564, -0.486, -0.402, -0.312, -0.215, -0.111, -0.000, NA, NA, NA, NA, NA, NA),
  c(-0.589, -0.523, -0.452, -0.374, -0.290, -0.200, -0.104, 0.000, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-0.418, -0.347, -0.270, -0.186, -0.097, 0.000, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-0.250, -0.173, -0.090, 0.000, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))

# Computes K1, K2 from omega1, omega2, using the lookup tables.
lineBinghamK1K2MLE <- function(omega1, omega2) {
  omega1s <- seq(from=0.02, to=0.32, by=0.02)
  omega2s <- seq(from=0.02, to=0.46, by=0.02)
  point <- matrix(c(omega1, omega2), 1, 2)
  k1 <- interp.surface(list(x=omega1s, y=omega2s, z=lineBinghamK1Table), point)
  k2 <- interp.surface(list(x=omega1s, y=omega2s, z=lineBinghamK2Table), point)
  c(k1, k2)
}

#' Maximum likelihood estimation of the Bingham distribution parameters.
#' 
#' Based on Tauxe (2010). May not work if the data are too concentrated or too dispersed? Results are great for n == 1000, so-so for n == 100, and poor for n == 10. The Bingham probability density is proportional to exp(-x^T A x), not exp(x^T A x). I prefer to normalize A so that tr A == 0.
#' @param us A list of lines.
#' @param weights A vector of real numbers (positive), of length equal to us. They need not sum to 1; the function automatically normalizes them to do so.
#' @return A list with members $a (symmetric 3x3 real matrix A), $values (a real 3D vector; the eigenvalues of A), and $vectors (a rotation matrix; the eigenvectors of A are the columns). The values are in descending order and sum to zero.
lineTauxeBinghamMLE <- function(us, weights=replicate(length(us), 1)) {
  # In eigen, the eigenvectors are descending and the eigenvectors are unit.
  #tBar <- arithmeticMean(lapply(us, function(u) outer(u, u)))
  w <- weights / sum(weights)
  ts <- lapply(1:length(us), function(i) {w[[i]] * outer(us[[i]], us[[i]])})
  tBar <- Reduce("+", ts)
  eig <- eigen(tBar, symmetric=TRUE)
  # Find MLE concentration parameters k1, k2, based on eigvals omega1 <= omega2 <= omega3.
  omega1 <- eig$values[[3]]
  omega2 <- eig$values[[2]]
  omega3 <- eig$values[[1]]
  k1k2 <- lineBinghamK1K2MLE(omega1, omega2)
  # Compute the eigenvalues of A. Normalize so that tr A == 0.
  values <- c(-k1k2[[1]], -k1k2[[2]], 0)
  values <- values - sum(values) / 3
  # Repackage into matrix A.
  vectors <- cbind(eig$vectors[,3], eig$vectors[,2], eig$vectors[,1])
  if (det(vectors) < 0)
    vectors[,3] <- -vectors[,3]
  a <- vectors %*% diag(values) %*% t(vectors)
  list(a=a, values=values, vectors=vectors)
}

# Test.
#vectors <- rotUniform()
#values <- c(8, 4, 0)
#a <- vectors %*% diag(values) %*% t(vectors)
#us <- lineBingham(a, n=1000)
#mle <- lineTauxeBinghamMLE(us)
#a
#mle$a
#vectors
#mle$vectors
#oriDistance(t(vectors), t(mle$vectors), group=oriLineInPlaneGroup)
#values - min(values)
#mle$values - min(mle$values)

#' 95% confidence region for the mean of the Bingham distribution.
#' 
#' This function sometimes fails, if the data set is too concentrated, dispersed, or small? This is not a particularly high-quality implementation of the technique.
#' @param ls A list of lines.
#' @param numPoints A real number (non-negative integer). If numPoints > 0, then this function constructs a curve for the boundary of the 95% confidence region.
#' @return A list with members $directions, $scatter, and $angles. The first two members are identical to $vectors and $values in lineProjectedMean. $angles is a pair of real numbers. They describe the 95% confidence region, as two distances from the mean toward the two other principal dispersions, measured in radians along the unit sphere's surface. If the inference fails, then the angles are NA. If numPoints > 0, then there is also a member $points, which is a list of numPoints + 1 lines delineating the confidence region.
lineBinghamInference <- function(us, numPoints=0) {
  # Compute eigensystem of (1 / n) SUM u_i u_i^T. Tauxe leaves out the n.
  n <- length(us)
  oriTsr <- arithmeticMean(lapply(us, function(u) outer(u, u)))
  eig <- eigen(oriTsr, symmetric=TRUE)
  # Find MLE concentration parameters k1, k2, based on eigvals omega1 <= omega2 <= omega3.
  omega1 <- eig$values[[3]]
  omega2 <- eig$values[[2]]
  omega3 <- eig$values[[1]]
  k1k2 <- lineBinghamK1K2MLE(omega1, omega2)
  # Use simple 95% confidence expressions from Tauxe Appendix C. She has another n here.
  k1 <- k1k2[[1]]
  k2 <- k1k2[[2]]
  epsilon32 <- 1.22 / (k2 * (omega2 - omega3))
  epsilon31 <- 1.22 / (k1 * (omega1 - omega3))
  result <- list(directions=eig$vectors,
       scatter=eig$values,
       angles=c(epsilon32, epsilon31))
  if (numPoints > 0) {
    f <- function(i) {
      theta <- 2 * pi * i / numPoints
      xyz <- c(1, epsilon32 * cos(theta), epsilon31 * sin(theta))
      rayNormalized(as.numeric(eig$vectors %*% xyz))
    }
    result$points <- lapply(0:numPoints, f)
  }
  result
}



### FITTING MIXTURES ###

#' Maximum likelihood estimation of a mixture of Bingham distributions.
#' 
#' Slow. Uses an iterative numerical optimization (expectation-minimization) wrapped around a numerical integration (to compute the normalizing constant) and sometimes another numerical optimization (to compute the single MLE) that contains a numerical integration (to compute the normalizing constant). Warning: Does not handle all of its internal errors carefully. More testing is warranted.
#' @param us A list of lines, of length n.
#' @param seeds A list of lines, of length k. Most importantly, k is the number of Bingham distributions to be mixed. No two seeds should be equal. Each seed should be the user's guess for where a Bingham mean might lie.
#' @param numSteps A real number (positive integer). A bound on the number of expectation-maximization iterations. Currently the algorithm always uses this many steps, so don't set it insanely high.
#' @return A list with members $phis, $as, $y, $assignments, $numSteps. $phis is a list of k real numbers (non-negative, summing to one), the weights on the distributions. $as is a list of real matrices (3x3, symmetric), the A parameters in the distributions. $numSteps is real (positive integer), the number of steps used. $assignments is a real vector (integers in [1, k]) indicating which distribution the us belong to, probabilistically. $y is a real matrix (k x n), a more detailed version of $assignments. The number yij is, roughly, the probability that uj belongs to distribution i.
lineMixedBinghamMLE <- function(us, seeds, numSteps, verbose=FALSE) {
  # Proclaim the distributions to be equally weighted.
  n <- length(us)
  k <- length(seeds)
  phis <- replicate(k, 1 / k)
  # Proclaim a distribution centered at each seed.
  f <- function(seed) {
    v <- rayOrthogonalUniform(seed)
    q <- cbind(seed, v, cross(seed, v))
    q %*% diag(c(-2, 1, 1)) %*% t(q)
  }
  as <- lapply(seeds, f)
  # Iteratively refine these arbitrary starting values for phis and as.
  step <- 0
  while (step < numSteps && 1 >= epsilon) {
    step <- step + 1
    # Compute the matrix y.
    cs <- sapply(as, lineBinghamIntegratedConstant)
    x <- sapply(1:n, function(j) sapply(1:k, function(i)
      phis[[i]] * cs[[i]] * exp(-us[[j]] %*% as[[i]] %*% us[[j]])))
    colsums <- sapply(1:n, function(j) sum(x[,j]))
    y <- sapply(1:n, function(j) sapply(1:k, function(i) 
      x[[i, j]] / colsums[[j]]))
    # Use y to update the weights phi and parameters A.
    phis <- sapply(1:k, function(i) mean(y[i,]))
    f <- function(i) {
      a <- lineTauxeBinghamMLE(us, weights=y[i,])$a
      if (Reduce("||", is.na(as.numeric(a))))
        a <- lineBinghamMLE(us, weights=y[i,], numSteps=10000)$a
      a
    }
    as <- lapply(1:k, f)
    if (verbose) {
      print(paste("step", as.character(step), "cs:"))
      print(cs)
      print(paste("step", as.character(step), "y:"))
      print(y)
      print(paste("step", as.character(step), "phis:"))
      print(phis)
      print(paste("step", as.character(step), "as:"))
      print(as)
    }
  }
  # Assign each line to a distribution probabilistically.
  ass <- sapply(1:n, function(j) which.max(as.numeric(rmultinom(1, 1, y[,j]))))
  list(phis=phis, as=as, y=y, assignments=ass, numSteps=step)
}

# Test. Because the above function sometimes fails, the test sometimes fails. Don't be surprised. For some reason the Tauxe MLE tends to work better than my MLE inside this mixed MLE, even though the opposite relationship is true outside the mixed MLE. And my integrated constant tends to work better than the saddlepoint-approximated constant again.
lineMixedBinghamMLETest <- function(ns, sd, numSteps) {
  # Make k matrices A.
  k <- length(ns)
  f <- function() {
    a <- rnorm(5, mean=0, sd=sd)
    rbind(c(a[[1]], a[[2]], a[[3]]), c(a[[2]], a[[4]], a[[5]]), c(a[[3]], a[[5]], -a[[1]] - a[[4]]))
  }
  as <- replicate(k, f(), simplify=FALSE)
  # Make synthetic data.
  uss <- lapply(1:k, function(i) lineBingham(as[[i]], ns[[i]]))
  us <- unlist(uss, recursive=FALSE)
  # Compute the mixed MLE, using admittedly fantastic seed means.
  mus <- lapply(uss, lineProjectedMean)
  print(date())
  mle <- lineMixedBinghamMLE(us, mus, numSteps=numSteps)
  print(date())
  # Make plots.
  lineEqualAreaPlot(us, colors=hues(as.numeric(sapply(1:k, function(i) replicate(ns[[i]], i)))))
  lineEqualAreaPlot(us, colors=hues(mle$assignments))
  list(as=as, mle=mle)
}
#test <- lineMixedBinghamMLETest(c(100, 100, 100), 4, 10)
#test$as[[1]]
#test$mle$as[[1]]



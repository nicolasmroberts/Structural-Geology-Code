


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### WATSON DISTRIBUTION ###

# Helper function for lineWatson. Returns one line sampled from the Watson distribution.
lineWatsonHelper <- function(f, bound, rot) {
  phi <- runif(1, min=0, max=pi)
  y <- runif(1, min=0, max=bound)
  while (y > f(phi)) {
    phi <- runif(1, min=0, max=pi)
    y <- runif(1, min=0, max=bound)
  }
  theta <- runif(1, min=-pi, max=pi)
  as.numeric(rot %*% cartesianFromSpherical(c(1, phi, theta)))
}

#' Normalizing constant for the Watson distribution.
#' 
#' That is, the number C such that C exp(kappa (mu^T u)^2) integrates to 1. By the way, C == 1 / 1F1(0.5, 1.5, kappa).
#' @param kappa A real number. The concentration parameter.
#' @return A real number.
lineWatsonNormalizer <- function(kappa) {
  if (kappa >= 0)
    (2 * sqrt(kappa)) / (sqrt(pi) * as.numeric(erfi(sqrt(kappa))))
  else
    (2 * sqrt(-kappa)) / (sqrt(pi) * as.numeric(erf(sqrt(-kappa))))
}

#' Sampling lines from the Watson distribution.
#' 
#' A naive acceptance-rejection sampling algorithm, based on bounding the density (with respect to the distance from mu) with a constant. For large kappa, this method grows inefficient. For kappa == 100, about 13 tries are needed per success. For kappa == -100, about 18 tries are needed.
#' @param mu A line. The mean of the distribution (if kappa > 0) or the pole to the girdle of the distribution (if kappa < 0).
#' @param kappa A real number. The concentration parameter.
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single line. If n is a positive integer, then a list of n lines.
lineWatson <- function(mu, kappa, n=NULL) {
  nu <- rayOrthogonalUniform(mu)
  rot <- cbind(cross(nu, mu), nu, mu)
  cc <- lineWatsonNormalizer(kappa)
  f <- function(phi) {exp(kappa * cos(phi)^2) * sin(phi) * cc / 2}
  if (kappa <= 0.5)
    bound <- cc / 2
  else
    bound <- exp(kappa - 1 / 2) * (2 * kappa)^(-1 / 2) * cc / 2
  if (is.null(n))
    lineWatsonHelper(f, bound, rot)
  else
    replicate(n, lineWatsonHelper(f, bound, rot), simplify=FALSE)
}

lineWatsonMLED3s <- c(0.001, seq(from=0.005, to=0.330, by=0.005), 0.333,
                      seq(from=0.34, to=0.99, by=0.01), 0.995, 0.999)
lineWatsonMLEKappaHats <- c(
  -500, -100, -50, -33.33, -25, -20, -16.67, -14.29, -12.25, -11.11, -9.992, -9.087, -8.327, -7.681, -7.126,
  -6.641, -6.215, -5.836, -5.495, -5.188, -4.908, -4.651, -4.415, -4.196, -3.993, -3.802, -3.624, -3.457,
  -3.298, -3.148, -3.006, -2.870, -2.741, -2.617, -2.499, -2.385, -2.275, -2.170, -2.068, -1.970, -1.874, -1.782, -1.692,
  -1.605, -1.520, -1.438, -1.357, -1.279, -1.202, -1.127, -1.053, -0.982, -0.911, -0.842, -0.774, -0.708,
  -0.642, -0.578, -0.514, -0.452, -0.390, -0.330, -0.270, -0.211, -0.152, -0.095, -0.038, -0.004,
  0.075, 0.184, 0.292, 0.398, 0.503, 0.606, 0.708, 0.809, 0.909, 1.008, 1.106, 1.204, 1.302, 1.399, 1.497,
  1.594, 1.692, 1.790, 1.888, 1.987, 2.087, 2.188, 2.289, 2.392, 2.496, 2.602, 2.709, 2.819,
  2.930, 3.044, 3.160, 3.280, 3.402, 3.529, 3.659, 3.764, 3.934, 4.079, 4.231, 4.389, 4.556, 4.731, 4.917,
  5.115, 5.326, 5.552, 5.797, 6.063, 6.354, 6.676, 7.035, 7.438, 7.897, 8.426, 9.043, 9.776,
  10.654, 11.746, 13.112, 14.878, 17.242, 20.560, 25.546, 33.866, 50.521, 100.510, 200.5, 1000.5)
lineWatsonMLEInterpolation <- approxfun(x=lineWatsonMLED3s, y=lineWatsonMLEKappaHats)

#' Maximum likelihood estimation of the Watson distribution parameters.
#' 
#' From Mardia and Jupp (2000, Section 10.3.2).
#' @param xs A list of lines.
#' @param shape NULL or character, either 'bipolar' or 'girdle'. If NULL, then this function chooses automatically.
#' @return A list with members $muHat (a line, the MLE of the mean), $kappaHat (a real number, the MLE of the concentration, $shape (character, either 'bipolar' or 'girdle'), $d3 (a positive real number, the D3 from which kappaHat was computed), and $eigenvalues (the eigenvalues of the T-bar matrix, in descending order).
lineWatsonMLE <- function(xs, shape=NULL) {
  # In eigen, the eigenvectors are descending and the eigenvectors are unit.
  tBar <- arithmeticMean(lapply(xs, function(x) outer(x, x)))
  eig <- eigen(tBar, symmetric=TRUE)
  # Pick a shape if necessary.
  if (is.null(shape)) {
    if (eig$values[[1]] - eig$values[[2]] >= eig$values[[2]] - eig$values[[3]])
      shape <- "bipolar"
    else
      shape <- "girdle"
  }
  # Shape determines which eigenvector is picked for muHat.
  if (shape == "bipolar")
    muHat <- eig$vectors[,1]
  else
    muHat <- eig$vectors[,3]
  # Get kappaHat from the lookup table.
  d3 <- as.numeric(muHat %*% tBar %*% muHat)
  if (d3 > 0.9)
    kappaHat <- 1 / (1 - d3)
  else if (d3 < 0.05)
    kappaHat <- -1 / (2 * d3)
  else
    kappaHat <- lineWatsonMLEInterpolation(d3)
  list(muHat=muHat, kappaHat=kappaHat, shape=shape, d3=d3, eigenvalues=eig$values)
}

# Test that the MLE and the sampling are both working. They work best when n is large. They work best when kappa is large, and somewhat better when kappa is positive than when kappa is negative.
lineWatsonTest <- function(kappa, n) {
  mu <- lineUniform()
  mle <- lineWatsonMLE(lineWatson(mu, kappa, n))
  c(mle$kappaHat, lineDistance(mu, mle$muHat))
}
#plot(t(replicate(100, lineWatsonTest(kappa=-10, n=10))), xlab="kappa", ylab="error in mean") # Often weird kappa > 0 results.
#plot(t(replicate(100, lineWatsonTest(kappa=-3, n=10))), xlab="kappa", ylab="error in mean") # Often weird kappa > 0.
#plot(t(replicate(100, lineWatsonTest(kappa=-1, n=10))), xlab="kappa", ylab="error in mean") # Often weird kappa > 0.
#plot(t(replicate(100, lineWatsonTest(kappa=1, n=10))), xlab="kappa", ylab="error in mean") # Often weird kappa < 0.
#plot(t(replicate(100, lineWatsonTest(kappa=3, n=10))), xlab="kappa", ylab="error in mean") # Often weird kappa < 0.
#plot(t(replicate(100, lineWatsonTest(kappa=10, n=10))), xlab="kappa", ylab="error in mean")
#plot(t(replicate(100, lineWatsonTest(kappa=-10, n=100))), xlab="kappa", ylab="error in mean")
#plot(t(replicate(100, lineWatsonTest(kappa=-3, n=100))), xlab="kappa", ylab="error in mean") # Often weird kappa > 0.
#plot(t(replicate(100, lineWatsonTest(kappa=-1, n=100))), xlab="kappa", ylab="error in mean") # Often weird kappa > 0.
#plot(t(replicate(100, lineWatsonTest(kappa=1, n=100))), xlab="kappa", ylab="error in mean") # Often weird kappa < 0.
#plot(t(replicate(100, lineWatsonTest(kappa=3, n=100))), xlab="kappa", ylab="error in mean")
#plot(t(replicate(100, lineWatsonTest(kappa=10, n=100))), xlab="kappa", ylab="error in mean")
#plot(t(replicate(100, lineWatsonTest(kappa=-10, n=1000))), xlab="kappa", ylab="error in mean")
#plot(t(replicate(100, lineWatsonTest(kappa=-3, n=1000))), xlab="kappa", ylab="error in mean")
#plot(t(replicate(100, lineWatsonTest(kappa=-1, n=1000))), xlab="kappa", ylab="error in mean")
#plot(t(replicate(100, lineWatsonTest(kappa=1, n=1000))), xlab="kappa", ylab="error in mean")
#plot(t(replicate(100, lineWatsonTest(kappa=3, n=1000))), xlab="kappa", ylab="error in mean")
#plot(t(replicate(100, lineWatsonTest(kappa=10, n=1000))), xlab="kappa", ylab="error in mean")

#' One-sample inference about the mean of the Watson distribution.
#' 
#' Assumes large concentration --- either kappa >> 0 or kappa << 0. From Mardia and Jupp (2000, Section 10.7.3).
#' @param xs A list of lines.
#' @param alpha A real number, between 0 and 1. The significance level for the confidence region.
#' @param shape NULL or character, either 'bipolar' or 'girdle'. If NULL, then this function chooses automatically.
#' @return A list with members $shape, $tBar, $rhs, $pvalue. $shape is either 'bipolar' or 'girdle'. If 'bipolar', then the confidence region consists of all lines u such that u^T %*% $tBar %*% u > $rhs. If 'girdle', then the confidence region consists of all lines u such that u^T %*% $tBar %*% u < $rhs. $pvalue is an R function that takes as input a line u0 and produces as output a real number in [0, 1] --- the p-value for the null hypothesis that the Watson mean is u0.
lineWatsonInference <- function(xs, alpha=0.05, shape=NULL) {
  # In eigen, the eigenvectors are descending and the eigenvectors are unit.
  tBar <- arithmeticMean(lapply(xs, function(x) outer(x, x)))
  eig <- eigen(tBar, symmetric=TRUE)
  # Pick a shape if necessary.
  if (is.null(shape)) {
    if (eig$values[[1]] - eig$values[[2]] >= eig$values[[2]] - eig$values[[3]])
      shape <- "bipolar"
    else
      shape <- "girdle"
  }
  n <-  length(xs)
  if (shape == "bipolar") {
    t1 <- eig$values[[1]]
    rhs <- t1 + (t1 - 1) * qf(alpha, 2, 2 * n - 2, lower.tail=FALSE) / (n - 1)
    func <- function(u0) {
      f <- as.numeric((t1 - u0 %*% tBar %*% u0) * (n - 1) / (1 - t1))
      1 - pf(f, 2, 2 * n - 2)
    }
  } else {
    t3 <- eig$values[[3]]
    rhs <- t3 * (1 + qf(alpha, 2, n - 2, lower.tail=FALSE) * 2 / (n - 2))
    func <- function(u0) {
      f <- as.numeric((u0 %*% tBar %*% u0 - t3) * (n - 2) / (2 * t3))
      1 - pf(f, 2, n - 2)
    }
  }
  list(shape=shape, tBar=tBar, rhs=rhs, pvalue=func)
}

# Tests. Works well for kappa >= 10 and n >= 10. Works well for kappa <= -10 and n >= 100. Works poorly in other cases, including all |kappa| == 1 cases.
lineWatsonInferenceExperiment <-  function(kappa, n, N) {
  f <- function(kappa, n) {
    mu <- lineUniform()
    us <- lineWatson(mu, kappa, n)
    inf <- lineWatsonInference(us)
    lhs <- as.numeric(mu %*% inf$tBar %*% mu)
    if (inf$shape == "bipolar")
      c(lhs - inf$rhs, inf$pvalue(mu))
    else
      c(inf$rhs - lhs, inf$pvalue(mu))
  }
  results <- replicate(N, f(kappa, n))
  coverageRateRegions <- sum(results[1,] > 0) / N
  coverageRatePvalues <- sum(results[2,] > 0.05) / N
  mismatches1 <- sum(results[1,] > 0 & results[2,] < 0.05)
  mismatches2 <- sum(results[1,] < 0 & results[2,] > 0.05)
  c(coverageRateRegions, coverageRatePvalues, mismatches1, mismatches2)
}
#lineWatsonInferenceExperiment(kappa=-100, n=1000, N=2000) # 0.943 0.943 0.000 0.000
#lineWatsonInferenceExperiment(kappa=-10, n=1000, N=2000) # 0.947 0.947 0.000 0.000
#lineWatsonInferenceExperiment(kappa=-1, n=1000, N=2000) # 0.7735 0.7735 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=1, n=1000, N=2000) # 0.8295 0.8295 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=10, n=1000, N=2000) # 0.969 0.969 0.000 0.000
#lineWatsonInferenceExperiment(kappa=100, n=1000, N=2000) # 0.9485 0.9485 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=-100, n=100, N=2000) # 0.9515 0.9515 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=-10, n=100, N=2000) # 0.9485 0.9485 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=-1, n=100, N=2000) # 0.5155 0.5155 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=1, n=100, N=2000) # 0.7095 0.7095 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=10, n=100, N=2000) # 0.96 0.96 0.00 0.00
#lineWatsonInferenceExperiment(kappa=100, n=100, N=2000) # 0.948 0.948 0.000 0.000
#lineWatsonInferenceExperiment(kappa=-100, n=30, N=2000) # 0.9085 0.9085 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=-10, n=30, N=2000) # 0.8795 0.8795 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=-1, n=30, N=2000) # 0.3645 0.3645 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=1, n=30, N=2000) # 0.5425 0.5425 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=10, n=30, N=2000) # 0.9565 0.9565 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=100, n=30, N=2000) # 0.962 0.962 0.000 0.000
#lineWatsonInferenceExperiment(kappa=-100, n=10, N=2000) # 0.636 0.636 0.000 0.000
#lineWatsonInferenceExperiment(kappa=-10, n=10, N=2000) # 0.5625 0.5625 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=-1, n=10, N=2000) # 0.4575 0.4575 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=1, n=10, N=2000) # 0.6255 0.6255 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=10, n=10, N=2000) # 0.9655 0.9655 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=100, n=10, N=2000) # 0.9465 0.9465 0.0000 0.0000

#' Multi-sample inference about the mean of the Watson distribution.
#' 
#' Assumes large sample sizes. Assumes that all of the data sets are Watson-distributed with the same unknown concentration kappa. From Mardia and Jupp (2000, Section 10.7.4).
#' @param xss A list of lists of lines.
#' @param shape NULL or character, either 'bipolar' or 'girdle'. If NULL, then this function chooses automatically.
#' @return A real number. The p-value for the null hypothesis that the means of the q distributions are all equal.
lineLargeMultiSampleWatsonInference <- function(xss, shape=NULL) {
  q <- length(xss)
  ns <- sapply(xss, length)
  tBars <- lapply(xss, function(xs) arithmeticMean(lapply(xs, function(x) outer(x, x))))
  n <- sum(ns)
  tBar <- arithmeticMean(lapply(unlist(xss, recursive=FALSE), function(x) outer(x, x)))
  eig <- eigen(tBar, symmetric=TRUE)
  # Pick a shape if necessary.
  if (is.null(shape)) {
    if (eig$values[[1]] - eig$values[[2]] >= eig$values[[2]] - eig$values[[3]])
      shape <- "bipolar"
    else
      shape <- "girdle"
  }
  if (shape == "bipolar") {
    t1 <- eig$vectors[,1]
    eHatT4 <- sum(sapply(xss, function(xs) sum(sapply(xs, function(x) dot(x, t1)^4)))) / n
    t1s <- sapply(tBars, function(mat) eigen(mat, symmetric=TRUE, only.values=TRUE)$values[[1]])
    t1 <- eig$values[[1]]
    chiSq <- (3 * t1 - 1) * (dot(ns, t1s) - n * t1) / (2 * (t1 - eHatT4))
    1 - pchisq(chiSq, 2 * q - 2)
  } else {
    t3 <- eig$vectors[,3]
    eHatT4 <- sum(sapply(xss, function(xs) sum(sapply(xs, function(x) dot(x, t3)^4)))) / n
    t3s <- sapply(tBars, function(mat) eigen(mat, symmetric=TRUE, only.values=TRUE)$values[[3]])
    t3 <- eig$values[[3]]
    chiSq <- (1 - 3 * t3) * (n * t3 - dot(ns, t3s)) / (2 * (t3 - eHatT4))
    1 - pchisq(chiSq, 2 * q - 2)
  }
}

#' Multi-sample inference about the mean of the Watson distribution.
#' 
#' Assumes tight concentration. Tests suggest that it works for |kappa| >= 10, but not for |kappa| == 1. Assumes that all of the data sets are Watson-distributed with the same unknown concentration kappa. From Mardia and Jupp (2000, Section 10.7.4).
#' @param xss A list of lists of lines.
#' @param shape NULL or character, either 'bipolar' or 'girdle'. If NULL, then this function chooses automatically.
#' @return A real number. The p-value for the null hypothesis that the means of the q distributions are all equal.
lineConcentratedMultiSampleWatsonInference <- function(xss, shape=NULL) {
  q <- length(xss)
  ns <- sapply(xss, length)
  tBars <- lapply(xss, function(xs) arithmeticMean(lapply(xs, function(x) outer(x, x))))
  n <- sum(ns)
  tBar <- arithmeticMean(lapply(unlist(xss, recursive=FALSE), function(x) outer(x, x)))
  eig <- eigen(tBar, symmetric=TRUE)
  # Pick a shape if necessary.
  if (is.null(shape)) {
    if (eig$values[[1]] - eig$values[[2]] >= eig$values[[2]] - eig$values[[3]])
      shape <- "bipolar"
    else
      shape <- "girdle"
  }
  if (shape == "bipolar") {
    t1s <- sapply(tBars, function(mat) eigen(mat, symmetric=TRUE, only.values=TRUE)$values[[1]])
    t1 <- eig$values[[1]]
    f <- (dot(ns, t1s) - n * t1) * (n - q) / (dot(ns, 1 - t1s) * (q - 1))
    1 - pf(f, 2 * q - 2, 2 * n - 2 * q)
  } else {
    t3s <- sapply(tBars, function(mat) eigen(mat, symmetric=TRUE, only.values=TRUE)$values[[3]])
    t3 <- eig$values[[3]]
    f <- (n * t3 - dot(ns, t3s)) * (n - 2 * q) / (dot(ns, t3s) * (2 * q - 2))
    1 - pf(f, 2 * q - 2, n - 2 * q)
  }
}

# Tests. The large-sample-size method is always too conservative, even at n = 3000? The high-concentration method works for |kappa| >= 10, for all sample sizes, but not for |kappa| == 1.
lineMultiSampleWatsonInferenceExperiment <-  function(kappa, n, N, q) {
  f <- function(kappa, n, q) {
    mu <- lineUniform()
    uss <- replicate(q, lineWatson(mu, kappa, n), simplify=FALSE)
    large <- lineLargeMultiSampleWatsonInference(uss)
    concen <- lineConcentratedMultiSampleWatsonInference(uss)
    c(large, concen)
  }
  results <- replicate(N, f(kappa, n, q))
  rateLarge <- sum(results[1,] > 0.05) / N
  rateConcen <- sum(results[2,] > 0.05) / N
  c(rateLarge, rateConcen)
}
#lineMultiSampleWatsonInferenceExperiment(kappa=-100, n=10, N=2000, q=2) # 0.9935 0.9490
#lineMultiSampleWatsonInferenceExperiment(kappa=-10, n=10, N=2000, q=2) # 0.9975 0.9460
#lineMultiSampleWatsonInferenceExperiment(kappa=-3, n=10, N=2000, q=2) # 0.9990 0.9335
#lineMultiSampleWatsonInferenceExperiment(kappa=-1, n=10, N=2000, q=2) # 0.999 0.866
#lineMultiSampleWatsonInferenceExperiment(kappa=1, n=10, N=2000, q=2) # 0.9995 0.8855
#lineMultiSampleWatsonInferenceExperiment(kappa=3, n=10, N=2000, q=2) # 0.9985 0.9615
#lineMultiSampleWatsonInferenceExperiment(kappa=10, n=10, N=2000, q=2) # 0.9995 0.9615
#lineMultiSampleWatsonInferenceExperiment(kappa=100, n=10, N=2000, q=2) # 0.997 0.953
#lineMultiSampleWatsonInferenceExperiment(kappa=-100, n=30, N=2000, q=2) # 0.9975 0.9515
#lineMultiSampleWatsonInferenceExperiment(kappa=-10, n=30, N=2000, q=2) # 0.9980 0.9505
#lineMultiSampleWatsonInferenceExperiment(kappa=-3, n=30, N=2000, q=2) # 0.9970 0.9285
#lineMultiSampleWatsonInferenceExperiment(kappa=-1, n=30, N=2000, q=2) # 0.9985 0.7645
#lineMultiSampleWatsonInferenceExperiment(kappa=1, n=30, N=2000, q=2) # 0.9985 0.8115
#lineMultiSampleWatsonInferenceExperiment(kappa=3, n=30, N=2000, q=2) # 0.9990 0.9605
#lineMultiSampleWatsonInferenceExperiment(kappa=10, n=30, N=2000, q=2) # 0.9965 0.9545
#lineMultiSampleWatsonInferenceExperiment(kappa=100, n=30, N=2000, q=2) # 0.9975 0.9460
#lineMultiSampleWatsonInferenceExperiment(kappa=-100, n=100, N=2000, q=2) # 0.999 0.961
#lineMultiSampleWatsonInferenceExperiment(kappa=-10, n=100, N=2000, q=2) # 0.9975 0.9480
#lineMultiSampleWatsonInferenceExperiment(kappa=-3, n=100, N=2000, q=2) # 0.9965 0.9220
#lineMultiSampleWatsonInferenceExperiment(kappa=-1, n=100, N=2000, q=2) # 0.9990 0.7505
#lineMultiSampleWatsonInferenceExperiment(kappa=1, n=100, N=2000, q=2) # 0.996 0.795
#lineMultiSampleWatsonInferenceExperiment(kappa=3, n=100, N=2000, q=2) # 0.9975 0.9565
#lineMultiSampleWatsonInferenceExperiment(kappa=10, n=100, N=2000, q=2) # 0.9985 0.9625
#lineMultiSampleWatsonInferenceExperiment(kappa=100, n=100, N=2000, q=2) # 0.9970 0.9495
#lineMultiSampleWatsonInferenceExperiment(kappa=-100, n=300, N=2000, q=2) # 0.9980 0.9485
#lineMultiSampleWatsonInferenceExperiment(kappa=-10, n=300, N=2000, q=2) # 0.9980 0.9555
#lineMultiSampleWatsonInferenceExperiment(kappa=-3, n=300, N=2000, q=2) # 0.9985 0.9410
#lineMultiSampleWatsonInferenceExperiment(kappa=-1, n=300, N=2000, q=2) # 0.9980 0.7735
#lineMultiSampleWatsonInferenceExperiment(kappa=1, n=300, N=2000, q=2) # 0.9975 0.8135
#lineMultiSampleWatsonInferenceExperiment(kappa=3, n=300, N=2000, q=2) # 0.9945 0.9655
#lineMultiSampleWatsonInferenceExperiment(kappa=10, n=300, N=2000, q=2) # 0.9975 0.9555
#lineMultiSampleWatsonInferenceExperiment(kappa=100, n=300, N=2000, q=2) # 0.9965 0.9550
#lineMultiSampleWatsonInferenceExperiment(kappa=-100, n=1000, N=2000, q=2) # 0.9975 0.9440
#lineMultiSampleWatsonInferenceExperiment(kappa=-10, n=1000, N=2000, q=2) # 0.9990 0.9575
#lineMultiSampleWatsonInferenceExperiment(kappa=-3, n=1000, N=2000, q=2) # 0.998 0.936
#lineMultiSampleWatsonInferenceExperiment(kappa=-1, n=1000, N=2000, q=2) # 0.9955 0.7835
#lineMultiSampleWatsonInferenceExperiment(kappa=1, n=1000, N=2000, q=2) # 0.9960 0.8115
#lineMultiSampleWatsonInferenceExperiment(kappa=3, n=1000, N=2000, q=2) # 0.9970 0.9625
#lineMultiSampleWatsonInferenceExperiment(kappa=10, n=1000, N=2000, q=2) # 0.9995 0.9625
#lineMultiSampleWatsonInferenceExperiment(kappa=100, n=1000, N=2000, q=2) # 0.9955 0.9530
#lineMultiSampleWatsonInferenceExperiment(kappa=-100, n=3000, N=2000, q=2) # 0.9990 0.9495 36minutes
#lineMultiSampleWatsonInferenceExperiment(kappa=-10, n=3000, N=2000, q=2) # 0.9995 0.9540
#lineMultiSampleWatsonInferenceExperiment(kappa=-3, n=3000, N=2000, q=2) # 0.9985 0.9260
#lineMultiSampleWatsonInferenceExperiment(kappa=-1, n=3000, N=2000, q=2) # 0.9975 0.7820
#lineMultiSampleWatsonInferenceExperiment(kappa=1, n=3000, N=2000, q=2) # 0.9965 0.8040
#lineMultiSampleWatsonInferenceExperiment(kappa=3, n=3000, N=2000, q=2) # 0.9975 0.9615
#lineMultiSampleWatsonInferenceExperiment(kappa=10, n=3000, N=2000, q=2) # 0.9965 0.9555
#lineMultiSampleWatsonInferenceExperiment(kappa=100, n=3000, N=2000, q=2) # 0.9990 0.9535



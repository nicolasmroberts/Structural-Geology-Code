


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# We should do more with the Fisher distribution:
# unbiased estimation of isotropic Fisher parameters (Mardia, 1972, p. 250)
# isotropic Fisher goodness of fit (Mardia, 1972, p. 252)
# uniformity of rays 256
# Watson-Williams test p. 214
# Watson-Williams two-sample test 263
# assumes Fisher concentrations equal, so also test that (266)



### FISHER DISTRIBUTION ###

# Generates one ray from the Fisher distribution centered at the third column of the rotation matrix, using accept-reject with f and its bound.
rayFisherHelper <- function(rot, f, bound) {
  # The azimuthal coordinate theta is uniform on the unit circle.
  theta <- runif(1, min=-pi, max=pi)
  # Perform acceptance-rejection sampling for phi.
  phi <- runif(1, min=0, max=pi)
  y <- runif(1, min=0, max=bound)
  while (y > f(phi)) {
    phi <- runif(1, min=0, max=pi)
    y <- runif(1, min=0, max=bound)
  }
  as.numeric(rot %*% cartesianFromSpherical(c(1, phi, theta)))
}

#' Sampling points from the Fisher distribution.
#' 
#' Also called von Mises-Fisher distribution or Langevin distribution. Uses expressions on p. 172 of Mardia and Jupp (2000), with a naive acceptance-rejection sampling algorithm. As kappa increases, this algorithm gets less and less efficient --- for example, about 19 tries per acceptance when kappa == 100.
#' @param mu A ray. The center of the distribution.
#' @param kappa A real number (positive). The concentration of the distribution.
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single ray. If n is a positive integer, then a list of n rays.
rayFisher <- function(mu, kappa, n=NULL) {
  # Make a rotation matrix with mu as its third column.
  nu <- rayOrthogonalUniform(mu)
  rot <- cbind(cross(nu, mu), nu, mu)
  # phi is distributed as f = ((k / (2 sinh k)) e^(k cos phi) sin phi) on [0, pi].
  f <- function(phi) {
    cosPhi <- cos(phi)
    if (kappa * cosPhi < 700)
      (kappa / (2 * sinh(kappa))) * exp(kappa * cos(phi)) * sin(phi)
    else
      kappa * sin(phi) / (exp(kappa * (1 - cosPhi)) - exp(kappa * (-1 - cosPhi)))
  }
  # Bound f.
  bb <- sqrt(sqrt(1 + 4 * kappa^2) - 1)
  if (kappa < 700) {
    aa <- exp((sqrt(1 + 4 * kappa^2) - 1) / 2)
    cc <- 2 * sqrt(2) * sinh(kappa)
    bound <- aa * bb / cc
  } else
    bound <- bb * exp(-1 / 2) / sqrt(2)
  if (is.null(n)) {
    rayFisherHelper(rot, f, bound)
  } else
    replicate(n, rayFisherHelper(rot, f, bound), simplify=FALSE)
}

# From Mardia and Jupp (2000), Appendix 3.2.
rayFisherMLEKappaHats <- c(
  0.000, 0.030, 0.060, 0.090, 0.120, 0.150, 0.180, 0.211, 0.241, 0.271,
  0.302, 0.332, 0.363, 0.394, 0.425, 0.456, 0.488, 0.519, 0.551, 0.583,
  0.615, 0.647, 0.680, 0.713, 0.746, 0.780, 0.814, 0.848, 0.883, 0.918,
  0.953, 0.989, 1.025, 1.062, 1.100, 1.137, 1.176, 1.215, 1.255, 1.295,
  1.336, 1.378, 1.421, 1.464, 1.508, 1.554, 1.600, 1.647, 1.696, 1.746,
  1.797, 1.849, 1.903, 1.958, 2.015, 2.074, 2.135, 2.198, 2.263, 2.330,
  2.401, 2.473, 2.549, 2.628, 2.711, 2.798, 2.888, 2.984, 3.085, 3.191,
  3.304, 3.423, 3.551, 3.687, 3.832, 3.989, 4.158, 4.341, 4.541, 4.759,
  4.998, 5.262, 5.555, 5.882, 6.250, 6.667, 7.143, 7.692, 8.333, 9.091,
  10.000, 11.111, 12.500, 14.286, 16.667, 20.000, 25.000, 33.333, 50.000, 100.000)
rayFisherMLEInterpolation <- approxfun(x=seq(from=0.00, to=0.99, by=0.01), y=rayFisherMLEKappaHats)

#' Maximum likelihood estimation of the Fisher parameters.
#' 
#' Based on Mardia and Jupp (2000, p. 198).
#' @param xs A list of rays.
#' @return A list with members $muHat (the mean ray, identical to rayProjectedMean), $rBar (non-negative real number), and $kappaHat (a positive real number).
rayFisherMLE <- function(xs) {
  xBar <- arithmeticMean(xs)
  rBar <- sqrt(dot(xBar, xBar))
  x0 <- xBar / rBar
  if (rBar > 0.9)
    kappaHat <- 1 / (1 - rBar)
  else
    kappaHat <- rayFisherMLEInterpolation(rBar)
  list(muHat=x0, rBar=rBar, kappaHat=kappaHat)
}

# Test.
#mu <- rayUniform()
#us <- rayFisher(mu=mu, kappa=10, n=100)
#rayFisherMLE(us)
#mu

# From Mardia and Jupp (2000, Appendix 3.1).
rayFisherConfidenceKappas <- c(
  seq(from=0.0, to=4.9, by=0.1),
  seq(from=5.0, to=9.8, by=0.2),
  seq(from=10.0, to=12.5, by=0.5),
  c(13.0, 14.0, 15.0, 20.0, 30.0, 40.0, 50.0, 100.0))
rayFisherConfidenceDeltas <- c(
  154.2, 152.9, 151.5, 150.0, 148.4, 146.6, 144.8, 142.8, 140.8, 138.6,
  136.3, 133.9, 131.4, 128.9, 126.2, 123.6, 120.9, 118.3, 115.6, 113.0,
  110.4, 107.9, 105.4, 103.1, 100.8, 98.6, 96.5, 94.5, 92.6, 90.8,
  89.0, 87.4, 85.8, 84.3, 82.8, 81.4, 80.1, 78.8, 77.6, 76.5,
  75.4, 74.3, 73.3, 72.3, 71.3, 70.4, 69.6, 68.7, 67.9, 67.1,
  66.4, 64.9, 63.6, 62.3, 61.1, 60.0, 58.9, 57.9, 56.9, 56.0,
  55.1, 54.3, 53.5, 52.7, 52.0, 51.3, 50.6, 50.0, 49.3, 48.7,
  48.2, 47.6, 47.1, 46.5, 46.0, 45.5, 44.4, 43.3, 42.3, 41.4,
  40.5, 39.7, 38.2, 36.8, 31.8, 25.8, 22.3, 19.9, 14.0) * degree
rayFisherConfidenceInterpolation <- approxfun(x=rayFisherConfidenceKappas, y=rayFisherConfidenceDeltas)

#' Confidence region for the Fisher distribution mean.
#' 
#' Theoretically better than rayFisherLargeSampleConfidence? But in my tests it seems no better. Coverage rates tend to be too small for small sample sizes n. n = 30 delivers around 93% coverage. n = 100 delivers adequate coverage. Based on Eq. (10.4.26) from Mardia and Jupp (2000).
#' @param xs A list of rays.
#' @return A list consisting of $angle (a real number in [0, pi]) and everything from rayFisherMLE.
rayFisherConfidence <- function(xs) {
  n <- length(xs)
  mle <- rayFisherMLE(xs)
  kappa <- n * mle$kappaHat * mle$rBar
  if (kappa > 100)
    mle$angle <- 140.2 * degree * kappa^(-0.5)
  else
    mle$angle <- rayFisherConfidenceInterpolation(kappa)
  mle
}

rayFisherConfidenceExperiment <- function(N, kappa, n, alpha=0.05) {
  f <- function(kappa, n, alpha) {
    mu <- rayUniform()
    us <- rayFisher(mu, kappa, n)
    uBar <- rayProjectedMean(us)
    delta <- rayFisherConfidence(us)$angle
    mu %*% uBar > cos(delta)
  }
  pHat <- sum(replicate(N, f(kappa, n, alpha))) / N
  se <- standardErrorProportion(N, pHat)
  c(pHat, pHat - 2 * se, pHat + 2 * se)
}
# These experiments show that kappa doesn't affect the accuracy much. But n = 30 is a bit small.
#rayFisherConfidenceExperiment(N=2000, kappa=1, n=10) # 0.8765000 0.8617862 0.8912138
#rayFisherConfidenceExperiment(N=2000, kappa=10, n=10) # 0.9005000 0.8871135 0.9138865
#rayFisherConfidenceExperiment(N=2000, kappa=100, n=10) # 0.9015000 0.8881735 0.9148265
#rayFisherConfidenceExperiment(N=2000, kappa=1, n=30) # 0.9270000 0.9153663 0.9386337
#rayFisherConfidenceExperiment(N=2000, kappa=10, n=30) # 0.9290000 0.9175144 0.9404856
#rayFisherConfidenceExperiment(N=2000, kappa=100, n=30) # 0.9240000 0.9121489 0.9358511
#rayFisherConfidenceExperiment(N=2000, kappa=1, n=100) # 0.9435000 0.9331745 0.9538255
#rayFisherConfidenceExperiment(N=2000, kappa=10, n=100) # 0.9405000 0.9299208 0.9510792
#rayFisherConfidenceExperiment(N=2000, kappa=100, n=100) # 0.9525000 0.9429875 0.9620125
#rayFisherConfidenceExperiment(N=2000, kappa=1, n=300) # 0.9490000 0.9391614 0.9588386
#rayFisherConfidenceExperiment(N=2000, kappa=10, n=300) # 0.9485000 0.9386159 0.9583841
#rayFisherConfidenceExperiment(N=2000, kappa=100, n=300) # 0.9525000 0.9429875 0.9620125

# This is for comparing our results to Allmendinger and Cardozo.
#mu <- rayUniform()
#u10s <- rayFisher(mu, 6, 10)
#geoTrendPlungeDegFileFromRays(u10s, "u10s.txt")
#u100s <- rayFisher(mu, 6, 100)
#geoTrendPlungeDegFileFromRays(u100s, "u100s.txt")
#rayFisherConfidence5(u10s)$angle / degree
#rayFisherLargeSampleConfidence(u10s)$angle / degree
#rayFisherConfidence(u100s)$angle / degree
#rayFisherLargeSampleConfidence(u100s)$angle / degree

#' Asymptotic confidence region for the Fisher distribution mean.
#' 
#' Assumes large sample size n. When the sample size is small, this function's reported angle is typically too small. n = 10 is too small (only 90% coverage). n = 30 is a bit small (about 94% coverage). Even n = 55 is borderline. But n = 70 is quite adequate. Uses Eq. (10.4.31) from Mardia and Jupp (2000).
#' @param xs A list of rays.
#' @param alpha A real number (in [0, 1]). The significance level --- for example, 0.05 for 95% confidence.
#' @return A list with members $muHat (a ray, identical to rayProjectedMean), $rBar (a non-negative real number), $angle (a real number in [0, pi]). $angle is the radius of the confidence region, in radians, measured along the surface of the sphere.
rayFisherLargeSampleConfidence <- function(xs, alpha=0.05) {
  xBar <- arithmeticMean(xs)
  rBar <- sqrt(dot(xBar, xBar))
  x0 <- xBar / rBar
  tMatrix <- arithmeticMean(lapply(xs, function(x) outer(x, x)))
  n <- length(xs)
  sinDelta <- sqrt(-log(alpha) * (1 - x0 %*% tMatrix %*% x0) / (n * rBar^2))
  list(muHat=x0, rBar=rBar, angle=arcSin(sinDelta))
}

rayFisherLargeSampleConfidenceExperiment <- function(N, kappa, n, alpha=0.05) {
  f <- function(kappa, n, alpha) {
    mu <- rayUniform()
    us <- rayFisher(mu, kappa, n)
    uBar <- rayProjectedMean(us)
    delta <- rayFisherLargeSampleConfidence(us, alpha)$angle
    mu %*% uBar > cos(delta)
  }
  pHat <- sum(replicate(N, f(kappa, n, alpha))) / N
  se <- standardErrorProportion(N, pHat)
  c(pHat, pHat - 2 * se, pHat + 2 * se)
}
# These experiments show that kappa doesn't matter much to the accuracy. But n = 30 is a bit small.
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=1, n=10) # 0.8945000 0.8807618 0.9082382
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=1, n=30) # 0.9385000 0.9277559 0.9492441
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=1, n=100) # 0.9530000 0.9435352 0.9624648
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=1, n=300) # 0.9510000 0.9413461 0.9606539
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=10, n=10) # 0.902500 0.889234 0.915766
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=10, n=30) # 0.9375000 0.9266747 0.9483253
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=10, n=100) # 0.9550000 0.9457291 0.9642709
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=10, n=300) # 0.9535000 0.9440832 0.9629168
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=100, n=10) # 0.9045000 0.8913562 0.9176438
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=100, n=30) # 0.9370000 0.9261344 0.9478656
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=100, n=100) # 0.9455000 0.9353482 0.9556518
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=100, n=300) # 0.9435000 0.9331745 0.9538255
# Even n = 55 may be a bit small.
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=1, n=55) # 0.9500000 0.9402532 0.9597468
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=10, n=55) # 0.9400000 0.9293793 0.9506207
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=100, n=55) # 0.9335000 0.9223575 0.9446425
# But n = 70 is quite adequate.
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=1, n=70) # 0.9465000 0.9364364 0.9565636
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=10, n=70) # 0.9485000 0.9386159 0.9583841
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=100, n=70) # 0.9475000 0.9375257 0.9574743

#' Confidence region for the Fisher distribution mean.
#' 
#' See Tauxe (2010, p. 214). Experiments with Fisher-distributed data sets suggest that the sample size n doesn't affect the accuracy much. kappa == 1 is too dispersed, but kappa == 3 is fine.
#' @param xs A list of rays.
#' @param alpha A real number (in [0, 1]). The significance level --- for example, 0.05 for 95% confidence.
#' @return A list with members $muHat (a ray, identical to rayProjectedMean), $kappaHat (a non-negative real number), $angle (a real number in [0, pi]). $angle is the radius of the confidence region, in radians, measured along the surface of the sphere.
rayFisherTauxe <- function(xs, alpha=0.05) {
  resultant <- Reduce("+", xs)
  r <- sqrt(dot(resultant, resultant))
  xBar <- resultant / r
  n <- length(xs)
  kappa <- (n - 1) / (n - r)
  angle <- arcCos(1 - (alpha^(1 / (1 - n)) - 1) * (n - r) / r)
  list(muHat=xBar, kappaHat=kappa, angle=angle)
}

rayFisherTauxeExperiment <- function(N, kappa, n, alpha=0.05) {
  f <- function(kappa, n, alpha) {
    mu <- rayUniform()
    us <- rayFisher(mu, kappa, n)
    tauxe <- rayFisherTauxe(us, alpha)
    mu %*% tauxe$muHat > cos(tauxe$angle)
  }
  pHat <- sum(replicate(N, f(kappa, n, alpha))) / N
  se <- standardErrorProportion(N, pHat)
  c(pHat, pHat - 2 * se, pHat + 2 * se)
}
# These experiments show that n doesn't affect the accuracy much. kappa == 1 is too dispersed, but kappa == 3 is fine.
#rayFisherTauxeExperiment(N=2000, kappa=1, n=5) # 0.9065000 0.8934802 0.9195198
#rayFisherTauxeExperiment(N=2000, kappa=3, n=5) # 0.9475000 0.9375257 0.9574743
#rayFisherTauxeExperiment(N=2000, kappa=10, n=5) # 0.9440000 0.9337176 0.9542824
#rayFisherTauxeExperiment(N=2000, kappa=30, n=5) # 0.9510000 0.9413461 0.9606539
#rayFisherTauxeExperiment(N=2000, kappa=100, n=5) # 0.9495000 0.9397072 0.9592928
#rayFisherTauxeExperiment(N=2000, kappa=1, n=10) # 0.9005000 0.8871135 0.9138865
#rayFisherTauxeExperiment(N=2000, kappa=3, n=10) # 0.9470000 0.9369809 0.9570191
#rayFisherTauxeExperiment(N=2000, kappa=10, n=10) # 0.9580000 0.9490294 0.9669706
#rayFisherTauxeExperiment(N=2000, kappa=30, n=10) # 0.9495000 0.9397072 0.9592928
#rayFisherTauxeExperiment(N=2000, kappa=100, n=10) # 0.9390000 0.9282968 0.9497032
#rayFisherTauxeExperiment(N=2000, kappa=1, n=30) # 0.8605000 0.8450055 0.8759945
#rayFisherTauxeExperiment(N=2000, kappa=3, n=30) # 0.9520000 0.9424401 0.9615599
#rayFisherTauxeExperiment(N=2000, kappa=10, n=30) # 0.9550000 0.9457291 0.9642709
#rayFisherTauxeExperiment(N=2000, kappa=30, n=30) # 0.9420000 0.9315467 0.9524533
#rayFisherTauxeExperiment(N=2000, kappa=100, n=30) # 0.9510000 0.9413461 0.9606539
#rayFisherTauxeExperiment(N=2000, kappa=1, n=100) # 0.8735000 0.8586341 0.8883659
#rayFisherTauxeExperiment(N=2000, kappa=3, n=100) # 0.9535000 0.9440832 0.9629168
#rayFisherTauxeExperiment(N=2000, kappa=10, n=100) # 0.9465000 0.9364364 0.9565636
#rayFisherTauxeExperiment(N=2000, kappa=30, n=100) # 0.9545000 0.9451802 0.9638198
#rayFisherTauxeExperiment(N=2000, kappa=100, n=100) # 0.9525000 0.9429875 0.9620125
#rayFisherTauxeExperiment(N=2000, kappa=1, n=300) # 0.8795000 0.8649412 0.8940588
#rayFisherTauxeExperiment(N=2000, kappa=3, n=300) # 0.9500000 0.9402532 0.9597468
#rayFisherTauxeExperiment(N=2000, kappa=10, n=300) # 0.9525000 0.9429875 0.9620125
#rayFisherTauxeExperiment(N=2000, kappa=30, n=300) # 0.9460000 0.9358922 0.9561078
#rayFisherTauxeExperiment(N=2000, kappa=100, n=300) # 0.9475000 0.9375257 0.9574743

# Tauxe's kappaHats are usually smaller than the other MLE's kappaHats (except when kappa == 1; then Tauxe's is often larger).
#us <- rayFisher(mu=rayUniform(), kappa=1, n=10)
#rayFisherMLE(us)$kappaHat
#rayFisherTauxe(us)$kappaHat

#' Invert the relationship between kappa-hat and alpha95.
#' 
#' Given n and alpha95, works back to the MLE kappa that produced that alpha95. Useful for processing results from the literature, where alpha95 is reported but not the (more fundamental) kappa. Assumes that the author used the same 'Tauxe' approximation as in rayFisherTauxe.
#' @param n A real number (positive integer). The sample size.
#' @param alpha95 A real number (positive). The alpha95 angle, in radians.
#' @return A real number (positive). The Fisher MLE kappa-hat.
rayFisherTauxeKappaHatFromAlpha95 <- function(n, alpha95) {
  ratio <- (1 - cos(alpha95)) / (1 - 0.05^(1 / (1 - n)))
  r <- n / (1 - ratio)
  kappaHat <- (n - 1) / (n - r)
  kappaHat
}

#' The relationship between alpha95 and kappa-hat.
#' 
#' Basically the inverse to rayFisherTauxeKappaHatFromAlpha95. Assumes the same kind of 'Tauxe' approximation.
#' @param n A real number (positive integer). The sample size.
#' @param kappaHat real number (positive). The Fisher MLE kappa-hat.
#' @return A real number (positive). The alpha95 angle, in radians.
rayFisherTauxeAlpha95FromKappaHat <- function(n, kappaHat) {
  r = n + (1 - n) / kappaHat
  alpha <- 0.05
  angle <- arcCos(1 - (alpha^(1 / (1 - n)) - 1) * (n - r) / r)
  angle
}



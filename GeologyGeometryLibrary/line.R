


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# A line is expressed as a unit 3D vector in Cartesian coordinates, as an R vector u = c(x, y, z) where x^2 + y^2 + z^2 == 1. It is crucial to remember that u and -u represent the same line. Many of the line functions here are implemented in terms of the ray functions of rays.R. In particular, when lines are tightly concentrated, you can often ignore their 'negative copies' and treat them as rays, with no appreciable effect on the statistics.



### MISCELLANY ###

#' The ray that represents the line and that is closest to center.
#' 
#' @param l A line (unit 3D vector).
#' @param center A ray (unit 3D vector).
lineNearestRepresentative <- function(l, center) {
  if (dot(l, center) < 0)
    -l
  else
    l
}

#' Distance between two lines, as the angle (in radians) between them.
#' 
#' @param u A line.
#' @param v A line.
#' @return A real number. The distance between the two lines, in radians. Between 0 and pi / 2, inclusive.
lineDistance <- function(u, v) {
  arcCos(abs(dot(u, v)))
}

#' L^2 variance of a set of lines about a point.
#' 
#' I'm not sure about the weighting on this. If you change it, change regression too.
#' @param us A list of lines.
#' @param center A line. Usually some kind of mean of the us.
#' @return A real number (in the interval [0, pi^2 / 8]).
lineVariance <- function(us, center) {
  sum(sapply(us, function(u) lineDistance(u, center)^2)) / (2 * length(us))
}

#' The Frechet (geodesic L^2) mean of a set of lines.
#'
#' An interative algorithm for computing the Frechet mean --- the line that minimizes the Frechet variance. The iterations continue until error squared of epsilon is achieved or numSteps iterations have been used. Try multiple seeds, to improve your chances of finding the global optimum.
#' @param us A list of lines.
#' @param numSeeds A real number (positive integer). How many us to try as seeds.
#' @param numSteps A real number (positive integer). Bound on how many iterations to use.
#' @return A list consisting of $variance (a real number), $mean (a line), $error (an integer) and $minEigenvalue (a real number). error should be 0 and minEigenvalue should be positive. Otherwise there was some problem in the optimization. If error is non-zero, then try increasing numSteps.
lineMeanVariance <- function(us, numSeeds=5, numSteps=100) {
  f <- function(phiTheta) {
    lineVariance(us, cartesianFromSpherical(c(1, phiTheta)))
  }
  seeds <- sample(us, numSeeds)
  # Line variance is never larger than pi^2 / 8. So this candidate will lose.
  best <- list(variance=(pi^2))
  for (seed in seeds) {
    # Find the minimum. Record it with diagnostics if it's the best.
    sol <- optim(sphericalFromCartesian(seed), f, hessian=TRUE, control=list(maxit=numSteps))
    if (sol$value < best$variance) {
      eigvals <- eigen(sol$hessian, symmetric=TRUE, only.values=TRUE)$values
      best <- list(variance=sol$value, mean=cartesianFromSpherical(c(1, sol$par)), error=sol$convergence, minEigenvalue=min(eigvals))
    }
  }
  best
}

#' Geodesic curve between two points on the unit sphere.
#' 
#' @param u A line.
#' @param v A line.
#' @param numSteps The number of line segments to be used in the approximation of the great circle arc from u to v.
#' @return A list of numSteps+1 lines. The first is u and the last is v. They are evenly spaced and in order.
lineGeodesicPoints <- function(u, v, numSteps=10) {
  if (dot(u, v) < 0)
    rayGeodesicPoints(u, -v, numSteps)
  else 
    rayGeodesicPoints(u, v, numSteps)
}



### EXTRINSIC METHODS ###

#' Projected mean and scatter of a set of lines.
#' 
#' @param us A list of lines.
#' @return A list with elements $values, $vectors. $vectors is a 3x3 real matrix whose columns are the principal directions of dispersion: mean, other direction along main girdle, pole to girdle. $values is a vector of three numbers, nonnegative, summing to 1, and decreasing. They capture the strength of the concentration about the principal directions.
lineMeanScatter <- function(us) {
  tMatrices <- lapply(us, function(u) outer(u, u))
  tMatrix <- arithmeticMean(tMatrices)
  eig <- eigen(tMatrix, symmetric=TRUE)
  eig
}

#' Shortcut convenience function for the projected mean.
#' 
#' @param us A list of lines.
#' @return A line.
lineProjectedMean <- function(us) {
  eig <- lineMeanScatter(us)
  eig$vectors[,1]
}

#' Bootstrapped projected mean with percentile confidence region and hypothesis tests.
#' 
#' The inference is based on percentiles of Mahalanobis distance in the tangent space at the mean of the bootstrapped means. The user should check that the bootstrapped means form a tight ellipsoidal cluster, before taking such a region seriously.
#' @param ls A list of lines.
#' @param numBoots A real number (positive integer). The number of bootstrapped means to compute. 10,000 might be a good value.
#' @param ... Other arguments to be passed to the underlying rayMahalanobisPercentiles function.
#' @return A list. See rayMahalanobisPercentiles for most of it. But the $pvalue function treats its input as a ray. An added $pvalueLine function treats its input as a line, so use that.
lineBootstrapInference <- function(ls, numBoots, ...) {
  boots <- replicate(numBoots, lineProjectedMean(sample(ls, length(ls), replace=TRUE)), simplify=FALSE)
  bootMean <- lineProjectedMean(boots)
  boots <- lapply(boots, function(u) {if (dot(u, bootMean) >= 0) u else -u})
  inf <- rayMahalanobisPercentiles(boots, bootMean, ...)
  inf$pvalueLine <- function(l) {
    if (dot(l, inf$center) < 0)
      inf$pvalue(-l)
    else
      inf$pvalue(l)
  }
  inf
}



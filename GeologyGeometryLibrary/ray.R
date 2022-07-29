


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# A ray is a unit 3D vector, expressed in Cartesian coordinates as an R vector of the form u = c(x, y, z), where x^2 + y^2 + z^2 == 1. A ray is tantamount to a line with a preferred direction along that line. This file offers R functions for dealing with rays. Some are preceded by detailed documentation. The others are not intended for use by end-users.



### MISCELLANY ###

#' Vector normalization (to have length 1).
#' 
#' @param v A d-dimensional vector. Cannot be of length 0.
#' @return A d-dimensional vector of length 1.
rayNormalized <- function(v) {
  v / sqrt(sum(v^2))
}

#' Projection of a vector onto the plane perpendicular to another vector.
#' 
#' For example, when v = c(0, 0, 1), returns the upward-most ray perpendicular to pole.
#' @param pole A d-dimensional vector perpendicular to the plane of interest. Must have non-zero length. Need not be unit.
#' @param v A d-dimensional vector. Should not be parallel to pole. Need not be unit.
#' @return A unit d-dimensional vector.
rayOrthogonalProjection <- function(pole, v) {
  rayNormalized(v - pole * dot(v, pole) / dot(pole, pole))
}

#' A ray perpendicular to a given vector, deterministically.
#' 
#' The result is deterministic but arbitrary. Due to the hairy ball theorem, the result cannot be chosen to depend smoothly on the input.
#' @param v A 3-dimensional vector. Need not be unit.
#' @return A ray, perpendicular to v.
rayOrthogonal <- function(v) {
  if (abs(v[1]) < 0.5)
    w <- c(1, 0, 0)
  else if (abs(v[2]) < 0.5)
    w <- c(0, 1, 0)
  else
    w <- c(0, 0, 1)
  rayOrthogonalProjection(v, w)
}

#' A ray perpendicular to a given vector, probabilistically.
#' 
#' In theory, the returned ray is uniformly chosen on the circle's worth of rays perpendicular to the given vector.
#' @param v A 3-dimensional vector. Need not be unit.
#' @return A ray, perpendicular to v.
rayOrthogonalUniform <- function(v) {
  rayOrthogonalProjection(v, rayUniform())
}

#' The great circle perpendicular to the given ray.
#' 
#' @param pole A ray.
#' @param numSteps The number of line segments to be used in the approximation of the great circle.
#' @return A list of numSteps+1 rays. The first and last one are identical. Otherwise, they are evenly spaced on the circle perpendicular to the given pole, and in order (either clockwise or counter-clockwise, depending on your viewpoint).
rayGreatCircle <- function(pole, numSteps=50) {
  v <- rayOrthogonal(pole)
  w <- cross(pole, v)
  angles <- (0:numSteps) * (2 * pi / numSteps)
  lapply(angles, function(a) {cos(a) * v + sin(a) * w})
}

#' A small circle with the given ray as its pole.
#' 
#' @param v A ray. The pole of the small circle.
#' @param r A real number. The radius of the circle around the pole, in radians, measured along the surface of the sphere. For example, the Arctic circle would be at radius r == 23.5 * degree from the north pole.
#' @return A list of numPoints+1 rays. The first and last one are identical. Otherwise, they are evenly spaced on the small circle, and in order (either clockwise or counter-clockwise, depending on your viewpoint).
raySmallCircle <- function(v, r, numPoints=50) {
  # Make a rotation matrix with v as its third column.
  v2 <- rayOrthogonal(v)
  rotation <- cbind(v2, cross(v, v2), v)
  # I'm not sure what the best way to handle NA values is.
  if (is.na(r))
    r <- 0
  # Make the small circle near the north pole and rotate it to near v.
  lapply(0:numPoints, function(s) as.numeric(rotation %*% cartesianFromSpherical(c(1, r, 2 * pi * s / numPoints))))
}

#' Geodesic between two points on the unit sphere.
#' 
#' Like rayGreatCircle, but specifying the great circle in terms of two points on it, rather than its pole.
#' @param u A ray.
#' @param v A ray.
#' @param numSteps The number of line segments to be used in the approximation of the great circle arc from u to v.
#' @return A list of numSteps+1 rays. The first is u and the last is v. They are evenly spaced and in order.
rayGeodesicPoints <- function(u, v, numSteps=10) {
  dotProduct <- dot(u, v)
  if (dotProduct >= 1)
    replicate(numSteps + 1, u, simplify=FALSE)
  else {
    angle <- arcCos(dotProduct)
    perp <- rayNormalized(v - dotProduct * u)
    lapply(0:numSteps, function(i) {cos(angle * i / numSteps) * u + sin(angle * i / numSteps) * perp})
  }
}

# Returns the ray that is fraction s of the way from u to v.
# Performs badly when u == -v.
rayInterpolation <- function(u, v, s) {
  dotProduct <- dot(u, v)
  if (dotProduct >= 1)
    u
  else {
    angle <- arcCos(dotProduct)
    perp <- rayNormalized(v - dotProduct * u)
    cos(s * angle) * u + sin(s * angle) * perp
  }
}

# Given a nonempty curve in the unit sphere (a list of unit 3D Cartesian vectors).
# Returns a list of curves, where each curve is entirely z >= 0 or entirely z <= 0.
# Whenever it crosses the z == 0 equator, ends and begins adjacent curves there.
rayCurvesUpperLower <- function(us) {
  curves <- list()
  signs <- c()
  curve <- list(us[[1]])
  currentSign <- sign(us[[1]][[3]])
  for (u in us[2:length(us)]) {
    if (sign(u[[3]]) == 0 || sign(u[[3]]) == currentSign)
      # Continue the current curve in the current hemisphere.
      curve[[length(curve) + 1]] <- u
    else if (currentSign == 0) {
      # Continue the current curve and commit to a hemisphere.
      curve[[length(curve) + 1]] <- u
      currentSign <- sign(u[[3]])
    } else {
      v <- curve[[length(curve)]]
      if (sign(v[[3]]) == 0) {
        # The current curve has already ended.
        curves[[length(curves) + 1]] <- curve
        signs[[length(signs) + 1]] <- currentSign
        # Start the new curve at v or -v followed by u.
        if (dot(v, u) > 0)
          curve <- list(v, u)
        else
          curve <- list(-v, u)
        currentSign <- sign(u[[3]])
      } else {
        # Find equatorial point between v and u.
        p <- rayNormalized(v + (v[[3]] / (v[[3]] - u[[3]])) * (u - v))
        # End current curve with whichever of p, -p is closer to v.
        if (dot(v, p) > 0)
          curve[[length(curve) + 1]] <- p
        else
          curve[[length(curve) + 1]] <- -p
        curves[[length(curves) + 1]] <- curve
        signs[[length(signs) + 1]] <- currentSign
        # Begin next curve with whichever is closer to u, followed by u.
        if (dot(u, p) > 0)
          curve <- list(p, u)
        else
          curve <- list(-p, u)
        currentSign <- sign(u[[3]])
      }
    }
  }
  # Finish the current curve and return.
  curves[[length(curves) + 1]] <- curve
  signs[[length(signs) + 1]] <- currentSign
  list(curves=curves, signs=signs)
}

# A list of four triangles that partition the unit sphere.
rayTetrahedron <- list(
  list(c(1, 1, 1) / sqrt(3), c(1, -1, -1) / sqrt(3), c(-1, 1, -1) / sqrt(3)),
  list(c(-1, -1, 1) / sqrt(3), c(1, -1, -1) / sqrt(3), c(-1, 1, -1) / sqrt(3)),
  list(c(-1, 1, -1) / sqrt(3), c(-1, -1, 1) / sqrt(3), c(1, 1, 1) / sqrt(3)),
  list(c(1, -1, -1) / sqrt(3), c(-1, -1, 1) / sqrt(3), c(1, 1, 1) / sqrt(3)))

# Given a triangle of unit vectors in counterclockwise order when viewed from outside.
# Returns a list of 4^numNonAdapt triangles of unit vectors, each in counterclockwise order again.
rayRefinedTriangle <- function(tri, numNonAdapt) {
  if (numNonAdapt == 0)
    list(tri)
  else {
    w1 <- rayNormalized(tri[[1]] + tri[[2]])
    w2 <- rayNormalized(tri[[2]] + tri[[3]])
    w3 <- rayNormalized(tri[[3]] + tri[[1]])
    unlist(list(
      rayRefinedTriangle(list(w1, w2, w3), numNonAdapt - 1),
      rayRefinedTriangle(list(tri[[1]], w1, w3), numNonAdapt - 1),
      rayRefinedTriangle(list(tri[[2]], w2, w1), numNonAdapt - 1),
      rayRefinedTriangle(list(tri[[3]], w3, w2), numNonAdapt - 1)),
      recursive=FALSE, use.names=FALSE)
  }
}

#' Triangular approximation to the unit sphere.
#' 
#' @param numNonAdapt A real number (non-negative integer). The number of times to refine the base triangulation of the sphere.
#' @return A list of triangles, where each triangle is a list of three rays. The number of triangles is 4^(1 + numNonAdapt), so each triangle has spherical area pi / 4^numNonAdapt.
rayTetrahedralSphere <- function(numNonAdapt) {
  unlist(
    lapply(rayTetrahedron, rayRefinedTriangle, numNonAdapt),
    recursive=FALSE, use.names=FALSE) 
}

#' Numerical integration of a scalar function on the sphere.
#' 
#' @param f An R function from {rays} to {real numbers}.
#' @param numNonAdapt A real number (non-negative integer). The number of times to refine the base triangulation of the sphere. Each increment of numNonAdapt increases the time and memory requirements by a factor of 4, but makes the approximation better.
#' @return The approximated integral of f over the sphere.
rayIntegral <- function(f, numNonAdapt=5) {
  tris <- rayTetrahedralSphere(numNonAdapt)
  grid <- lapply(tris, function(tri) rayNormalized(tri[[1]] + tri[[2]] + tri[[3]]))
  area <- pi / 4^numNonAdapt
  sum(sapply(grid, f)) * area
}

# Tests.
#rayIntegral(function(u) 1, 5) # Should be 4 * pi. Yep.
#rayIntegral(function(u) u[[3]], 5) # Notice u[[3]] == cos(phi). Should be 0. Yep.
#rayIntegral(function(u) abs(u[[3]]), 5) # Should be 2 * pi. Numerical integration underestimates.
#rayIntegral(function(u) exp(-u[[3]]), 5) # Should be 2 * pi * (e - 1 / e). Yep.

#' Uniformly random points on the unit sphere.
#' 
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single ray. If n is a positive integer, then a list of n rays.
rayUniform <- function(n=NULL) {
  if (is.null(n))
    cartesianFromHorizontal(c(runif(1, 0, 2 * pi), runif(1, -1, 1)))
  else
    replicate(n, rayUniform(), simplify=FALSE)
}

#' Fitting a small circle to some points on the unit sphere.
#' 
#' @param us A list of rays.
#' @return A list consisting of $pole (a ray, the pole to the small circle), $angle (the distance from the pole to the small circle, in radians), $error (0 if and only if minimization succeeds; if 1, then increase numSteps), and $minEigenvalue (worry about the result if this isn't positive).
rayBestFitSmallCircle <- function(us, numSeeds=5, numSteps=1000) {
  f <- function(phiTheta) {
    p <- cartesianFromSpherical(c(1, phiTheta))
    angles <- sapply(us, function(u) arcCos(dot(p, u)))
    var(angles)
  }
  # Find the minimum, starting from a few seeds.
  best <- list(value=(2 * pi^2))
  for (i in 1:numSeeds) {
    seed <- sphericalFromCartesian(rayUniform())[2:3]
    solution <- optim(seed, f, lower=c(0, -pi), upper=c(pi, pi), hessian=TRUE, control=list(maxit=numSteps), method="L-BFGS-B")
    if (solution$value <= best$value)
      best <- solution
  }
  # Report diagnostic information.
  eigvals <- eigen(best$hessian, symmetric=TRUE, only.values=TRUE)$values
  p <- cartesianFromSpherical(c(1, best$par))
  angles <- sapply(us, function(u) arcCos(dot(p, u)))
  a <- mean(angles)
  if (a > pi / 2) {
    a <- pi - a
    p <- -p
  }
  list(pole=p, angle=a, error=best$convergence, minEigenvalue=min(eigvals))
}

# Test with exact or inexact data. Works well.
#pole <- rayUniform()
#angle <- runif(1, 0, pi)
#us <- raySmallCircle(pole, angle, numPoints=10)
##us <- lapply(us, function(u) rayNormalized(u + rnorm(3, 0, 0.1)))
#rayFit <- rayBestFitSmallCircle(us, numSeeds=10)
#rayFit
#rayEqualAreaPlot(us, curves=list(raySmallCircle(rayFit$pole, rayFit$angle)))



### EXTRINSIC METHODS ###

#' Extrinsic mean and scalar scatter of a set of rays.
#' 
#' Scatter varies between 0 (tight concentration) and 1 (wide dispersion). This scatter is denoted 1 - R-bar in Mardia and Jupp (2000), p. 163. Arguably the preferred measure of scatter should be 2 * (1 - R-bar) or 1 - R-bar^2, but this function implements neither of those.
#' @param us A list of rays.
#' @return A list with elements $mean (ray) and $scatter (a real number between 0 and 1). 
rayMeanScatter <- function(us) {
  u <- arithmeticMean(us)
  r <- sqrt(dot(u, u))
  list(mean=(u / r), scatter=(1 - r))
}

#' Extrinsic mean of a set of rays.
#' 
#' Convenience shortcut for rayMeanScatter(us)$mean.
#' @param us A list of rays.
#' @return A ray.
rayProjectedMean <- function(us) {
  rayNormalized(arithmeticMean(us))
}

#' Bootstrapped extrinsic mean with percentile confidence region and hypothesis tests.
#' 
#' The inference is based on percentiles of Mahalanobis distance in the tangent space at the mean of the bootstrapped means. The user should check that the bootstrapped means form a tight ellipsoidal cluster, before taking such a region seriously.
#' @param ls A list of rays.
#' @param numBoots A real number (positive integer). The number of bootstrapped means to compute. 10,000 might be a good value.
#' @param ... Other arguments to be passed to the underlying rayMahalanobisPercentiles function.
#' @return A list. See rayMahalanobisPercentiles.
rayBootstrapInference <- function(ls, numBoots, ...) {
  boots <- replicate(numBoots, rayProjectedMean(sample(ls, length(ls), replace=TRUE)), simplify=FALSE)
  bootMean <- rayProjectedMean(boots)
  rayMahalanobisPercentiles(boots, bootMean, ...)
}



### DIFFERENTIAL GEOMETRY ###

#' Distance between two points on the unit sphere.
#' 
#' @param u A ray.
#' @param v A ray.
#' @return A real number. The angular difference between the two rays, in radians. Always between 0 and pi, inclusive.
rayDistance <- function(u, v) {
  arcCos(dot(u, v))
}

#' L^2 variance of a set of rays about a point.
#' 
#' I'm not sure about the weighting on this. If you change it, change regression too.
#' @param us A list of rays.
#' @param center A ray. Usually some kind of mean of the us.
#' @return A real number (in the interval [0, pi^2 / 2]).
rayVariance <- function(us, center) {
  sum(sapply(us, function(u) rayDistance(u, center)^2)) / (2 * length(us))
}

#' The Frechet (geodesic L^2) mean of a set of rays.
#'
#' An interative algorithm for computing the Frechet mean --- the ray that minimizes the Frechet variance. The iterations continue until error squared of epsilon is achieved or numSteps iterations have been used. Try multiple seeds, to improve your chances of finding the global optimum.
#' @param us A list of rays.
#' @param numSeeds A real number (positive integer). How many us to try as seeds.
#' @param numSteps A real number (positive integer). Bound on how many iterations to use.
#' @return A list consisting of $variance (a real number), $mean (a ray), $error (an integer) and $minEigenvalue (a real number). error should be 0 and minEigenvalue should be positive. Otherwise there was some problem in the optimization. If error is non-zero, then try increasing numSteps.
rayMeanVariance <- function(us, numSeeds=5, numSteps=100) {
  f <- function(phiTheta) {
    rayVariance(us, cartesianFromSpherical(c(1, phiTheta)))
  }
  seeds <- sample(us, numSeeds)
  # Ray variance is never larger than pi^2 / 2. So this candidate will lose.
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

#' The exponential map on the unit sphere.
#' 
#' Regards v as a vector in the tangent space of the unit sphere at the point p. Returns a point q on the unit sphere, a distance of |v| from p, in the direction of v. Partial inverse to rayLog.
#' @param p A ray.
#' @param v A 3-dimensional vector, perpendicular to p, not necessarily unit.
#' @return A ray.
rayExp <- function(p, v) {
  normV <- sqrt(dot(v, v))
  if (normV == 0)
    p
  else {
    r <- rbind(p, v / normV, cross(p, v / normV))
    as.numeric(t(r) %*% c(cos(normV), sin(normV), 0))
  }
}

#' Inverse exponential map on the unit sphere.
#' 
#' Returns a vector v in the tangent space to the unit sphere at p, such that rayExp(p, v) = q.
#' @param p A ray.
#' @param q A ray.
#' @return A 3-dimensional vector, perpendicular to p, not necessarily unit.
rayLog <- function(p, q) {
  normV <- arcCos(dot(p, q))
  w <- rayNormalized(cross(p, q))
  v <- normV * cross(w, p)
  if (dot(v, q) >= 0)
    v
  else
    -v
}

# Tests that those are inverses.
#p <- rayUniform()
#v <- runif(1, min=0.1, max=3.1) * rayOrthogonalUniform(p)
#sqrt(dot(v, v))
#rayDistance(p, rayExp(p, v))
#v
#rayLog(p, rayExp(p, v))

#' Wrapping a plane around the unit sphere.
#' 
#' This function composes the exponential map with a non-canonical isomorphism to the plane R^2. The isomorphism is specified by the user through a rotation matrix R. The first row of R is regarded as a point p on the unit sphere. The other two rows define an isomorphism between R^2 and the tangent plane to the unit sphere at p. w is mapped through this isomorphism into the tangent plane, and then into the sphere via the exponential map. Inverse to rayTangentVectorFromPoint.
#' @param w A 2-dimensional real vector.
#' @param rotation A 3x3 real matrix (special orthogonal).
#' @return A ray.
rayPointFromTangentVector <- function(w, rotation) {
  v <- as.numeric(t(rotation) %*% c(0, w))
  q <- rayExp(rotation[1,], v)
  q
  # Does this other code do the same thing?
  #v <- c(sqrt(1 - vec[[1]]^2 - vec[[2]]^2), vec)
  #as.numeric(t(rotation) %*% v)
}

#' Unwrapping unit sphere into a tangent plane.
#' 
#' Inverse to rayPointFromTangentVector.
#' @param q A ray.
#' @param rotation A 3x3 real matrix (special orthogonal).
#' @return A 2-dimensional real vector.
rayTangentVectorFromPoint <- function(q, rotation) {
  v <- rayLog(rotation[1,], q)
  w <- as.numeric(rotation %*% v)[2:3]
  w
  # Does this faster code do the same thing?
  #as.numeric(rotation %*% q)[2:3]
}

# These tests demonstrate that the preceding two functions are inverses.
#p <- rayUniform()
#q <- rayUniform()
#perp <- rayOrthogonalUniform(p)
#rot <- rbind(p, perp, cross(p, perp))
#q
#w <- rayTangentVectorFromPoint(q, rot)
#rayPointFromTangentVector(w, rot)
#w
#rayTangentVectorFromPoint(rayPointFromTangentVector(w, rot), rot)

#' Principal component analysis in the tangent space.
#'
#' Appropriate only if the sample is tightly concentrated near the center.
#' @param us A list of rays.
#' @param center A ray. Typically the Frechet mean of the us.
#' @param numPoints A real number (integer, 0 or >= 3). The number of points to return on each of the four geodesics through the center.
#' @return A list consisting of $rotation (3x3 rotation matrix), $magnitudes (2D real vector, nonnegative), $directions (2x2 real matrix, whose columns are unit-length vectors), $pcsFromRay (an R function from rays to 2-dimensional vectors), and $rayFromPCs (an R function from 2-dimensional vectors to rays). The $magnitudes are in decreasing order. The $directions are the corresponding directions, suitable for use in rayPointFromTangentVector along with $rotation. If numPoints >= 1, then there is also a $curves field (list of two lists of (2 numPoints + 1) rays).
rayPrincipalComponentAnalysis <- function(us, center, numPoints=0) {
  # Map the rays into the tangent space at the mean.
  perp <- rayOrthogonalUniform(center)
  rot <- rbind(center, perp, cross(center, perp))
  vs <- lapply(us, rayTangentVectorFromPoint, rot)
  # Compute the usual covariance stuff.
  covar <- arithmeticMean(lapply(vs, function(v) {outer(v, v)}))
  eig <- eigen(covar, symmetric=TRUE)
  mags <- sqrt(eig$values)
  dirs <- eig$vectors
  # Compute functions for transferring back and forth.
  pcsFromRay <- function(u) {
    as.numeric(t(dirs) %*% rayTangentVectorFromPoint(u, rot))
  }
  rayFromPCs <- function(w) {
    rayPointFromTangentVector(as.numeric(dirs %*% w), rot)
  }
  result <- list(rotation=rot, magnitudes=mags, directions=dirs, pcsFromRay=pcsFromRay, rayFromPCs=rayFromPCs)
  # Include points for visualization, if desired.
  if (numPoints >= 1) {
    f <- function(s, magDir) {
      rayPointFromTangentVector(s / numPoints * magDir, rot)
    }
    curve1 <- lapply(-numPoints:numPoints, f, mags[[1]] * dirs[,1])
    curve2 <- lapply(-numPoints:numPoints, f, mags[[2]] * dirs[,2])
    result$curves <- list(curve1, curve2)
  }
  result
}

#' Elliptical region based on percentiles of Mahalanobis distance in the tangent space.
#' 
#' @param us A list of rays.
#' @param center A ray. Usually something like the Frechet mean of us.
#' @param alpha A real number, which is assumed to be 0.05 unless you specify 0.01.
#' @param numPoints A real number (non-negative integer). The resolution with which to sample the boundary curve of the (1 - alpha) * 100% percentile region.
#' @param doIsotropic Logical. If TRUE, forces the inverse covariance to the identity matrix and hence the region to be circular.
#' @return A list. $us is us. $center is center. $covarInv is the inverse covariance matrix in the tangent space, which is just the identity if doIsotropic is TRUE. rotation is the rotation used in rayTangentVectorFromPoint, etc. $q000, $q025, $q050, $q075, $q095, $q099 are quantiles of Mahalanobis norm. $pvalue is an R function that assigns to any ray its p-value, meaning the fraction of us that are farther from center than the given ray. $alpha is alpha. $angles is a pair of two angles, in radians, giving the semi-axis lengths of the confidence ellipse. If numPoints > 0, then the list also has an element $points. $points is a list of numPoints + 1 rays describing the boundary of the region, in order, with the first and last points identical.
rayMahalanobisPercentiles <- function(us, center, alpha=0.05, numPoints=0, doIsotropic=FALSE) {
  # Map the rays into the tangent space at the mean.
  perp <- rayOrthogonalUniform(center)
  rot <- rbind(center, perp, cross(center, perp))
  vs <- lapply(us, rayTangentVectorFromPoint, rot)
  # Compute the usual covariance stuff.
  if (doIsotropic)
    covarInv <- diag(c(1, 1))
  else {
    covar <- arithmeticMean(lapply(vs, function(v) {outer(v, v)}))
    covarInv <- solve(covar)
  }
  norms <- sapply(vs, function(v) {sqrt(v %*% covarInv %*% v)})
  empiricalCDF <- ecdf(norms)
  # Build the p-value function.
  f <- function(u) {
    v <- rayTangentVectorFromPoint(u, rot)
    1 - empiricalCDF(sqrt(v %*% covarInv %*% v))
  }
  # Compute some basic information about the confidence region.
  qs <- quantile(norms, probs=c(0.00, 0.25, 0.50, 0.75, 0.95, 0.99, 1.00), names=FALSE)
  if (alpha != 0.05 && alpha != 0.01)
    alpha <- 0.05
  if (alpha == 0.01)
    q <- qs[[6]]
  else
    q <- qs[[5]]
  eig <- eigen(covarInv, symmetric=TRUE)
  radii <- q * eig$values^(-1 / 2)
  # Record most of the results.
  result <- list(
    us=us, pvalue=f, center=center, covarInv=covarInv, rotation=rot,
    q000=qs[[1]], q025=qs[[2]], q050=qs[[3]], q075=qs[[4]], q095=qs[[5]], 
    q099=qs[[6]], q100=qs[[7]], alpha=alpha, angles=radii)
  if (numPoints > 0) {
    # Make the ellipse in the 2D tangent space.
    circle <- lapply(
      0:numPoints, 
      function(s) c(cos(s * 2 * pi / numPoints), sin(s * 2 * pi / numPoints)))
    vs <- lapply(circle, function(v) as.numeric(eig$vectors %*% (radii * v)))
    # Embed the tangent space in 3D and wrap it.
    result$points <- lapply(vs, rayPointFromTangentVector, rot)
    result$alpha <- alpha
  }
  result
}



### KENT DISTRIBUTION ###

# There is some code in kent.R, but some of it works quite poorly, so I probably haven't released it.



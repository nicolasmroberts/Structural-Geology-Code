


# Copyright 2016-2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# The complete 3D orientation of an object can be described using three perpendicular unit vectors in right-handed order. And actually the third vector is redundant, as it is simply the cross product of the other two. For example, the orientation of a duck in flight could be described by one vector pointing from its heart toward its head, and a second vector pointing from its heart out of its back. The particular convention that you use doesn't matter, as long as you use it consistently. You put these three vectors into the rows of a matrix R, which is then special orthogonal. This matrix is the rotation that rotates the vectors back to the x-, y-, and z-coordinate axes. In this file, 'rotation matrix' means 'real 3x3 matrix that is special orthogonal'.



### SYMMETRY GROUPS ###

# Frequently orientations are subject to some symmetry group G, which is a finite set of rotations (satisfying certain properties). For any rotation Q in G, the rotation matrix Q %*% R represents the same orientation as R does. Here are some common symmetry groups.

#' When trivial symmetry is used, the orientations are simply rotations. This case is so important that we have separate code, in rotations.R, for doing it. But let's include the trivial group, for completeness.
oriTrivialGroup <- list(diag(c(1, 1, 1)))

#' Ray-in-plane symmetry is applicable to faults-with-slip-directions, such as slickensides (and certain minerals).
oriRayInPlaneGroup <- list(
  diag(c(1, 1, 1)),
  diag(c(-1, 1, -1)))

#' Line-in-plane symmetry is applicable to foliation-lineations, cylindrical fold orientations, triaxial ellipsoid orientations, and earthquake focal mechanisms (and olivine).
oriLineInPlaneGroup <- list(
  diag(c(1, 1, 1)),
  diag(c(1, -1, -1)),
  diag(c(-1, 1, -1)),
  diag(c(-1, -1, 1)))

#' Trigonal trapezohedral is the point group of alpha-quartz.
oriTrigonalTrapezohedralGroup <- list(
  diag(c(1, 1, 1)),
  rotMatrixAboutZ(pi * 2 / 3),
  rotMatrixAboutZ(pi * 4 / 3), 
  diag(c(1, -1, -1)),
  rotMatrixAboutZ(pi * 2 / 3) %*% diag(c(1, -1, -1)) %*% t(rotMatrixAboutZ(pi * 2 / 3)),
  rotMatrixAboutZ(pi * 4 / 3) %*% diag(c(1, -1, -1)) %*% t(rotMatrixAboutZ(pi * 4 / 3)))

#' Hexagonal trapezohedral is the point group of beta-quartz.
oriHexagonalTrapezohedralGroup <- list(
  diag(c(1, 1, 1)),
  rotMatrixAboutZ(pi * 1 / 3),
  rotMatrixAboutZ(pi * 2 / 3),
  diag(c(-1, -1, 1)), 
  rotMatrixAboutZ(pi * 4 / 3), 
  rotMatrixAboutZ(pi * 5 / 3), 
  diag(c(1, -1, -1)),
  rotMatrixAboutZ(pi * 1 / 3) %*% diag(c(1, -1, -1)) %*% t(rotMatrixAboutZ(pi * 1 / 3)),
  rotMatrixAboutZ(pi * 2 / 3) %*% diag(c(1, -1, -1)) %*% t(rotMatrixAboutZ(pi * 2 / 3)),
  diag(c(-1, -1, 1)) %*% diag(c(1, -1, -1)) %*% diag(c(-1, -1, 1)),
  rotMatrixAboutZ(pi * 4 / 3) %*% diag(c(1, -1, -1)) %*% t(rotMatrixAboutZ(pi * 4 / 3)),
  rotMatrixAboutZ(pi * 5 / 3) %*% diag(c(1, -1, -1)) %*% t(rotMatrixAboutZ(pi * 5 / 3)))

#' Replicating a list of orientations into their multiple rotation representatives.
#' 
#' @param rs A list of rotation matrices. The representative rotations Rs.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A list of rotation matrices. The set G Rs. For each orientation G R, its |G| representatives are spread out in this list.
oriSymmetrizedRotations <- function(rs, group) {
  unlist(lapply(group, function(g) lapply(rs, function(r) g %*% r)), recursive=FALSE, use.names=FALSE)
}



### MISCELLANEOUS METHODS ###

#' The distance between two orientations as points in SO(3) / G.
#'
#' @param r A rotation matrix.
#' @param q A rotation matrix.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A real number (in the interval [0, pi]). The distance from G R to G R.
oriDistance <- function(r, q, group) {
  qRT <- tcrossprod(q, r)
  trGQRT <- max(sapply(group, function(g) tr(g %*% qRT)))
	arcCos((trGQRT - 1) / 2)
}

#' Diameter of a set of orientations.
#' 
#' @param rs A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A real number (non-negative).
oriDiameter <- function(rs, group) {
  f <- function(i) {
    max(sapply(1:(i - 1), function(j) oriDistance(rs[[i]], rs[[j]], group)))
  }
  max(sapply(2:(length(rs)), f))
}

#' Selecting an orientation representative near a given rotation.
#' 
#' @param r A rotation matrix.
#' @param center A rotation matrix.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A rotation matrix. The element of G R that is closest to center.
oriNearestRepresentative <- function(r, center, group) {
  gRs <- lapply(group, function(g) g %*% r)
  trCTGRs <- sapply(gRs, function(gr) tr(crossprod(center, gr)))
  i <- which.max(trCTGRs)
  gRs[[i]]
}

#' Selecting orientation representatives near a given rotation.
#' 
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A list of rotation matrices. For each R, the element of G R that is closest to center.
oriNearestRepresentatives <- function(rs, center=rs[[1]], group) {
  lapply(rs, oriNearestRepresentative, center, group)
}

#' The Frechet (geodesic L^2) variance of a set of orientations.
#'
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix. Often the mean of the rs.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A real number (non-negative). The variance of the points G R about the point G center.
oriVariance <- function(rs, center, group) {
  dists <- sapply(rs, oriDistance, center, group)
  sum(dists^2) / (2 * length(rs))
}

#' The Frechet (geodesic L^2) mean of a set of orientations as points in SO(3) / G.
#'
#' An iterative algorithm for computing the Frechet mean --- the orientation that minimizes the Frechet variance. The iterations continue until change-squared of epsilon is achieved or numSteps iterations have been used.
#' @param rs A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param numSeeds A real number (positive integer). How many seeds to try, in trying to find a global optimum.
#' @param numSteps A real number (positive integer). Bound on how many iterations to use.
#' @return A list consisting of $mean (a special orthogonal real 3x3 matrix), $variance (a real number), $changeSquared (a real number), and $numSteps (a non-negative integer). changeSquared is the square of the size of the final step. numSteps is the number of iterations used.
oriMeanVariance <- function(rs, group, numSeeds=5, numSteps=1000) {
  seeds <- sample(rs, numSeeds)
  # No variance is ever larger than pi^2 / 2 < 5.
  best <- list(5)
  for (seed in seeds) {
    rBar <- seed
    changeSquared <- epsilon + 1.0
    k <- 0
    while (changeSquared >= epsilon && k < numSteps) {
      w <- diag(c(0, 0, 0))
      for (r in rs) {
        rBarTGRs <- lapply(group, function(g) crossprod(rBar, g %*% r))
        i <- which.max(sapply(rBarTGRs, tr))
        w <- w + rotLog(rBarTGRs[[i]])
      }
      w <- w / length(rs)
      rBar <- rBar %*% rotExp(w)
      changeSquared <- tr(crossprod(w, w))
      k <- k + 1
    }
    var <- oriVariance(rs, rBar, group)
    if (var < best[[1]])
      best <- list(var, rBar, changeSquared, k)
  }
  list(variance=best[[1]], mean=best[[2]], changeSquared=best[[3]], numSteps=best[[4]])
}

#' The Frechet (geodesic L^2) mean. Convenience shortcut for oriMeanVariance.
#' 
#' @param rs A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param ... Additional parameters to be passed to oriMeanVariance, such as numSeeds and numSteps.
#' @return A rotation matrix. The Frechet mean.
oriFrechetMean <- function(rs, group, ...) {
  oriMeanVariance(rs, group, ...)$mean
}



### BOOTSTRAPPING ###

#' Inference for mean orientations, based on non-parametric bootstrapping.
#' 
#' This function bootstraps the Frechet orientation mean, returning various pieces of information. The first is a list of the bootstrapped means. The user should oriEqualAnglePlot and oriEqualVolumePlot this list, to make sure that the means form a fairly tight ellipsoidal cluster. If so, then the second piece of information may be used: An R function that takes an orientation R0 as input, as produces as output a p-value for the hypothesis that the mean of the population is R0. This p-value is the fraction of means that are farther from the mean of means than R0 is, based on the Mahalanobis distance of the means.
#' @param rs A list of rotation matrices. The data.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param ... Additional parameters to be passed to the underlying oriFrechetMean, such as numSeeds and numSteps.
#' @return A list consisting of $bootstraps (a list of rotation matrices), $pvalueOri (an R function from {orientation matrices} to {real numbers}), and everything returned by rotMahalanobisPercentiles.
oriBootstrapInference <- function(rs, numBoots, group, ...) {
  # This lapply used to be a replicate, but replicate works crazily with ....
  print("fraction complete:")
  f <- function(i) {
    print(i / numBoots)
    oriFrechetMean(rs=sample(rs, length(rs), replace=TRUE), group=group, ...)
  }
  boots <- lapply(1:numBoots, f)
  bootMean <- oriFrechetMean(boots, group)
  boots <- oriNearestRepresentatives(boots, bootMean, group)
  infer <- rotMahalanobisPercentiles(boots, bootMean)
  pfunc <- function(r) {
    max(sapply(group, function(g) infer$pvalue(g %*% r)))
  }
  infer$pvalueOri <- pfunc
  infer$bootstraps <- boots
  infer
}



### TANGENT SPACE METHODS ###

# The tangent space methods all operate in the tangent space at a point. To apply them to orientations, simply pick the representative for each orientation that is closest to the point, using oriNearestRepresentative(s). As always, tangent space methods are best for data that are tightly concentrated.

# When working modulo a symmetry group, there is an additional upper bound on how spread out the data can be. Let r be the minimum distance between any two group elements. In other words, r is the distance from the identity to the nearest other group element. See rotSeparation. Let d be the diameter of the data set (oriDiameter). Then we can safely work with representatives as long as d < r / 2. For then oriNearestRepresentative never leaves the chosen set of representatives.



### DEPRECATED ###

#' Projected arithmetic mean of a set of orientations as points in SO(3) / G.
#'
#' An iterative algorithm for computing the projected arithmetic mean of orientations rather than rotations. See Bachmann et al. (2010). The iterations continue until the change is small or the allowed number of iterations has been exhausted.
#' @param rs A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param numSeeds A real number (positive integer). The number of seeds to try, in trying to find a global optimum.
#' @param numSteps A real number (positive integer). Bound on how many iterations to use.
#' @return A list consisting of $values (as in rotMeanScatter), $rotations (as in rotMeanScatter), $changeSquared (a real number), and $numSteps (a non-negative integer). changeSquared is the square of the size of the final step. numSteps is the number of iterations used.
oriMeanScatter <- function(rs, group, numSeeds=1, numSteps=1000) {
  seeds <- sample(rs, numSeeds)
  best <- c(-2)
  for (seed in seeds) {
    rBar <- seed
    # No concentration is ever bigger than 1.
    concenOld <- -6
    concenNew <- -4
    k <- 0
    while ((concenNew - concenOld)^2 >= epsilon && k < numSteps) {
      rsNear <- oriNearestRepresentatives(rs, rBar, group)
      meanScatter <- rotMeanScatter(rsNear)
      concenOld <- concenNew
      concenNew <- meanScatter$values[[1]]
      rBar <- meanScatter$rotations[[1]]
      k <- k + 1
    }
    if (concenNew > best[[1]])
      best <- list(concenNew, meanScatter, (concenNew - concenOld)^2, k)
  }
  list(values=best[[2]]$values, rotations=best[[2]]$rotations, changeSquared=best[[3]], numSteps=best[[4]])
}

#' Projected mean orientation. Convenience shortcut for oriMeanScatter.
#' 
#' @param rs A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param numSeeds A real number (positive integer). The number of seeds to try, in trying to find a global optimum.
#' @param numSteps A real number (positive integer). Bound on how many iterations to use.
#' @return A rotation matrix.
oriProjectedMean <- function(rs, group, numSeeds=1, numSteps=1000) {
  oriMeanScatter(rs, group, numSeeds, numSteps)$rotations[[1]]
}



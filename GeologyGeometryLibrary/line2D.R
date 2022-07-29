


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# Some of the machinery for 3D lines still works: lineMeanScatter, lineProjectedMean, lineDistance. Even lineWellner and lineWellnerInference.



### LINES IN 2D ###

#' Uniformly random 2D lines.
#' 
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single 2D line. If n is a positive integer, then a list of n 2D lines.
line2DUniform <- ray2DUniform

#' Random 2D lines, drawn from the wrapped normal distribution on the unit circle.
#' 
#' @param mean A 2D line.
#' @param sd A real number (positive). The standard deviation sigma of the underlying normal distribution, in radians.
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single 2D line. If n is a positive integer, then a list of n 2D lines.
line2DWrappedNormal <- ray2DWrappedNormal

#' Bootstrapped extrinsic mean for 2D lines.
#' 
#' Essentially ray2DBootstrapInference, adapted for lines.
#' @param ls A list of 2D lines.
#' @param numBoots A real number (positive integer). The number of bootstrapped means to compute. 10,000 might be a good value.
#' @return A list ($center, $bootstraps, $pvalue, $q000, $q025, $q050, $q075, $q095, $q099, $q100). $bootstraps are the bootstraps. $center is their mean. The other fields are quantiles of distance from the mean, in radians, among the bootstraps. For example, a 95% confidence region consists of all 2D rays within $q095 of $center. $pvalue is an R function from 2D rays to real numbers, assigning a p-value to any given null hypothesis for the mean.
line2DBootstrapInference <- function(ls, numBoots) {
  boots <- replicate(numBoots, lineProjectedMean(sample(ls, length(ls), replace=TRUE)), simplify=FALSE)
  bootMean <- lineProjectedMean(boots)
  dists <- sapply(boots, lineDistance, bootMean)
  empiricalCDF <- ecdf(dists)
  # Build the p-value function.
  f <- function(u) {
    1 - empiricalCDF(lineDistance(u, bootMean))
  }
  # Compute a few popular percentiles.
  qs <- quantile(dists, probs=c(0.00, 0.25, 0.50, 0.75, 0.95, 0.99, 1.00), names=FALSE)
  list(center=bootMean, bootstraps=boots, pvalue=f,
       q000=qs[[1]], q025=qs[[2]], q050=qs[[3]], q075=qs[[4]], q095=qs[[5]], q099=qs[[6]], q100=qs[[7]])
}

#' Rose plot for 2D lines.
#' 
#' Nearly identical to ray2DRosePlot, but with each line represented by two rays. In theory this rose diagram should be perfectly symmetric about the origin. In practice, data points often fall on the boundary between bins and are resolved arbitrarily, producing some leakage from bins into neighboring bins and hence asymmetry?
#' @param lines A list of 2D or 3D real vectors. In each one, only the 2D projection (the first two components) is used. It must be non-zero but need not be unit.
#' @param weights A vector of real numbers. The weights to attach to the lines.
#' @param ... Other options to be passed to the underlying rayRosePlot.
#' @return NULL.
line2DRosePlot <- function(lines, weights=replicate(length(lines), 1), ...) {
  rays <- c(lines, lapply(lines, function(v) -v))
  wghts <- c(weights / 2, weights / 2)
  ray2DRosePlot(rays, wghts, ...)
}



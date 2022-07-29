


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### RAYS IN 2D ###

# Some of the general machinery still works: rayNormalized, rayMeanScatter, rayProjectedMean, rayDistance.

#' Uniformly random points on the unit circle.
#' 
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single 2D ray. If n is a positive integer, then a list of n 2D rays.
ray2DUniform <- function(n=NULL) {
  if (is.null(n)) {
    angle <- runif(1, min=-pi, max=pi)
    c(cos(angle), sin(angle))
  }
  else
    replicate(n, rayUniform(), simplify=FALSE)
}

#' Random points on the unit circle, drawn from the wrapped normal distribution.
#' 
#' @param mean A 2D ray.
#' @param sd A real number (positive). The standard deviation sigma of the underlying normal distribution, in radians.
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single 2D ray. If n is a positive integer, then a list of n 2D rays.
ray2DWrappedNormal <- function(mean, sd, n=NULL) {
  if (is.null(n))
    ray2DWrappedNormal(mean, sd, 1)[[1]]
  else {
    angle <- atan2(mean[[2]], mean[[1]])
    angles <- rnorm(n, mean=angle, sd=sd)
    lapply(angles, function(a) c(cos(a), sin(a)))
  }
}

#' Bootstrapped extrinsic mean for 2D rays.
#' 
#' The inference is based on percentiles of distance from the mean of the bootstrapped means. The user should check that the bootstrapped means form a unimodal symmetric distribution, before taking such a region seriously.
#' @param ls A list of 2D rays.
#' @param numBoots A real number (positive integer). The number of bootstrapped means to compute. 10,000 might be a good value.
#' @return A list ($center, $bootstraps, $pvalue, $q000, $q025, $q050, $q075, $q095, $q099, $q100). $bootstraps are the bootstraps. $center is their mean. The other fields are quantiles of distance from the mean, in radians, among the bootstraps. For example, a 95% confidence region consists of all 2D rays within $q095 of $center. $pvalue is an R function from 2D rays to real numbers, assigning a p-value to any given null hypothesis for the mean.
ray2DBootstrapInference <- function(ls, numBoots) {
  boots <- replicate(numBoots, rayProjectedMean(sample(ls, length(ls), replace=TRUE)), simplify=FALSE)
  bootMean <- rayProjectedMean(boots)
  dists <- sapply(boots, rayDistance, bootMean)
  empiricalCDF <- ecdf(dists)
  # Build the p-value function.
  f <- function(u) {
    1 - empiricalCDF(rayDistance(u, bootMean))
  }
  # Compute a few popular percentiles.
  qs <- quantile(dists, probs=c(0.00, 0.25, 0.50, 0.75, 0.95, 0.99, 1.00), names=FALSE)
  list(center=bootMean, bootstraps=boots, pvalue=f,
       q000=qs[[1]], q025=qs[[2]], q050=qs[[3]], q075=qs[[4]], q095=qs[[5]], q099=qs[[6]], q100=qs[[7]])
}

#' Rose plot for 2D rays.
#' 
#' The data are binned. (Binning starts at the positive x-axis and proceeds counterclockwise. In a future release maybe we could add an offset, to start the binning elsewhere.) Each datum can be given a weight, which is tantamount to repeating the datum in the data set. Each petal in the Rose plot represents the total weight in a bin, either by area or by length. Concentric circles indicate the 10%, 20%, 30%, etc. weight levels.
#' @param rays A list of 2D or 3D real vectors. In each one, only the 2D projection (the first two components) is used. It must be non-zero but need not be unit.
#' @param weights A vector of real numbers. The weights to attach to the rays.
#' @param numBins The number of bins to use. For example, numBins == 36 means 10-degree-wide bins.
#' @param inner A real number (non-negative). This parameter reserves a circle of empty space at the center of the plot. More precisely, the plot takes place between radius inner and radius 1. inner == 0 seems traditional in geology, but I feel that inner == 0.25 (say) is prettier and easier to read.
#' @param areal Logical. If FALSE, then each petal's length is proportional to its bin's weight. If TRUE, then each petal's area is proportional to its bin's weight. areal == FALSE seems traditional in geology, but all data visualization advice I've ever seen suggests that areal == TRUE is the right choice.
#' @return NULL.
ray2DRosePlot <- function(rays, weights=replicate(length(rays), 1), numBins=36, inner=0, areal=TRUE) {
  # Bin the data, filling each bin with its fraction of the total weight.
  bins <- replicate(numBins, 0)
  for (i in 1:length(rays)) {
    heading <- atan2(rays[[i]][[2]], rays[[i]][[1]]) %% (2 * pi)
    bin <- heading %/% (2 * pi / numBins) + 1
    bins[[bin]] <- bins[[bin]] + weights[[i]]
  }
  bins <- bins / sum(weights)
  # Prepare to scale bins and benchmarks.
  if (areal)
    radius <- function(w) {sqrt(w * (1 - inner^2) + inner^2)}
  else
    radius <- function(w) {w * (1 - inner) + inner}
  # Draw benchmarks.
  plot.new()
  plot.window(xlim=c(-1, 1), ylim=c(-1, 1))
  xs <- cos((0:360) * (2 * pi / 360))
  ys <- sin((0:360) * (2 * pi / 360))
  for (w in seq(from=0, to=1, by=0.1))
    lines(radius(w) * xs, radius(w) * ys)
  # Draw the petals.
  for (bin in 1:numBins) {
    headings <- seq(from=((bin - 1) * 2 * pi / numBins), to=(bin * 2 * pi / numBins), by=(2 * pi / 360))
    xs <- radius(bins[[bin]]) * cos(headings)
    ys <- radius(bins[[bin]]) * sin(headings)
    heading <- bin * 2 * pi / numBins
    r <- radius(max(bins[[bin]], bins[[(bin %% numBins) + 1]]))
    xs <- c(xs, inner * cos(heading), r * cos(heading))
    ys <- c(ys, inner * sin(heading), r * sin(heading))
    lines(xs, ys)
  }
}



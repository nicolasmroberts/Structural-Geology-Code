


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### ORIENTATION PLOTS ###

#' X-Z-X Euler angle plot of orientations as symmetric sets of rotations.
#'
#' This function is much like rotEulerAnglePlot, but it plots |G| symmetric copies of the points, curves, and triangles. This plot does not have an equal-angle or equal-volume property. Indeed, it is extremely distorted near certain gimbal lock lines. Curves and triangles are not supported.
#' @param points A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param showBoundary Logical. Whether to show the bounding box.
#' @param ... Plotting options to be passed to plot3D. If colors is used, its length should equal that of points.
#' @return NULL.
oriEulerAnglePlot <- function(points, group, showBoundary=FALSE, ...) {
  rotEulerAnglePlot(points=oriSymmetrizedRotations(points, group), showBoundary=showBoundary, ...)
}

#' Axis-angle plot of orientations as symmetric sets of rotations.
#' 
#' This function is much like rotAxisAnglePlot, but it plots |G| symmetric copies of the points, curves, and triangles. Curves that cross the boundary are automatically clipped. On the other hand, triangles are not clipped. Thus you probably don't want to use any triangles.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param ... Plotting options to be passed to rotBallPlot. If colors is used, its length should equal that of points.
#' @return NULL.
oriAxisAnglePlot <- function(group, points=list(), curves=list(), triangles=list(), ...) {
  pointss <- oriSymmetrizedRotations(points, group)
  f <- function(g, ps) {lapply(ps, function(p) g %*% p)}
  curvess <- Reduce(c, lapply(group, function(g) lapply(curves, function(ps) f(g, ps))))
  triangless <- Reduce(c, lapply(group, function(g) lapply(triangles, function(ps) f(g, ps))))
  rotAxisAnglePlot(pointss, curvess, triangless, ...)
}

#' Equal-angle plot of orientations as symmetric sets of rotations.
#'
#' This function is much like rotEqualAnglePlot, but it plots |G| symmetric copies of the points, curves, and triangles.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param ... Plotting options to be passed to rotBallPlot. If colors is used, its length should equal that of points.
#' @return NULL.
oriEqualAnglePlot <- function(group, points=list(), curves=list(), triangles=list(), ...) {
  pointss <- oriSymmetrizedRotations(points, group)
  f <- function(g, ps) {lapply(ps, function(p) g %*% p)}
  curvess <- Reduce(c, lapply(group, function(g) lapply(curves, function(ps) f(g, ps))))
  triangless <- Reduce(c, lapply(group, function(g) lapply(triangles, function(ps) f(g, ps))))
  rotEqualAnglePlot(pointss, curvess, triangless, ...)
}

#' Equal-volume plot of orientations as symmetric sets of rotations.
#' 
#' This function is much like rotEqualVolumePlot, but it plots |G| symmetric copies of the points, curves, and triangles. Curves that cross the boundary are automatically clipped. On the other hand, triangles are not clipped. Thus you probably don't want to use any triangles.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param ... Plotting options to be passed to rotBallPlot. If colors is used, its length should equal that of points.
#' @return NULL.
oriEqualVolumePlot <- function(group, points=list(), curves=list(), triangles=list(), ...) {
  pointss <- oriSymmetrizedRotations(points, group)
  f <- function(g, ps) {lapply(ps, function(p) g %*% p)}
  curvess <- Reduce(c, lapply(group, function(g) lapply(curves, function(ps) f(g, ps))))
  triangless <- Reduce(c, lapply(group, function(g) lapply(triangles, function(ps) f(g, ps))))
  rotEqualVolumePlot(pointss, curvess, triangless, ...)
}

#' Rodrigues plot of orientations as symmetric sets of rotations.
#'
#' This function is much like rotRodriguesPlot, but it plots |G| symmetric copies of the points, curves, and triangles.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param ... Plotting options to be passed to rotBallPlot. If colors is used, its length should equal that of points.
#' @return NULL.
oriRodriguesPlot <- function(group, points=list(), curves=list(), triangles=list(), ...) {
  pointss <- oriSymmetrizedRotations(points, group)
  f <- function(g, ps) {lapply(ps, function(p) g %*% p)}
  curvess <- Reduce(c, lapply(group, function(g) lapply(curves, function(ps) f(g, ps))))
  triangless <- Reduce(c, lapply(group, function(g) lapply(triangles, function(ps) f(g, ps))))
  rotRodriguesPlot(pointss, curvess, triangless, ...)
}

#' Equal-angle plot of geodesic regression results.
#' 
#' The user may choose to draw extra spurs, joining the Rs to their corresponding predictions on the curve.
#' @param xs A vector of real numbers.
#' @param rs A list of rotation matrices.
#' @param regr The output of oriGeodesicRegression with those xs and rs.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param spurs A list of rotation matrices. When not NULL, it's usually equal to Rs.
#' @param ... Other plotting options to be passed to rotEqualAnglePlot.
#' @return NULL.
oriGeodesicRegressionPlot <- function(xs, rs, regr, group, colors=hues(xs), spurs=NULL, ...) {
  rotA <- rotExp(min(xs) * regr$m) %*% regr$b
  rotB <- rotExp(mean(range(xs)) * regr$m) %*% regr$b
  rotC <- rotExp(max(xs) * regr$m) %*% regr$b
  curves <- list(rotGeodesicPoints(rotA, rotB, 20), rotGeodesicPoints(rotB, rotC, 20))
  if (!is.null(spurs)) {
    spurCurves <- lapply(1:length(xs), function(i) {rotGeodesicPoints(spurs[[i]], regr$prediction(xs[[i]]), 10)})
    curves <- c(curves, spurCurves)
  }
  oriEqualAnglePlot(group=group, points=rs, curves=curves, colors=colors, ...)
}



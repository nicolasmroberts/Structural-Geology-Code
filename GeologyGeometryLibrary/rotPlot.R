


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### BASIC PLOTTING ###

#' x-z-x Euler angle plot of the space of rotations.
#'
#' The plot represents each rotation as a point in Euler angle space, which is a (2 pi) x (pi) x (2 pi) box. This plot does not have an equal-angle or equal-volume property. Indeed, it is extremely distorted near certain gimbal lock lines. Curves and triangles are not supported.
#' @param points A list of rotation matrices.
#' @param showBoundary Logical. Whether to show the bounding box.
#' @param ... Additional plotting options to be passed to plot3D.
#' @return NULL.
rotEulerAnglePlot <- function(points, showBoundary=FALSE, ...) {
  xzxs <- lapply(points, rotXZXAnglesFromMatrix)
  if (showBoundary)
    curves <- list(
      list(c(pi, pi, pi), c(-pi, pi, pi), c(-pi, 0, pi), c(pi, 0, pi), c(pi, pi, pi)),
      list(c(pi, pi, -pi), c(-pi, pi, -pi), c(-pi, 0, -pi), c(pi, 0, -pi), c(pi, pi, -pi)),
      list(c(pi, pi, pi), c(pi, pi, -pi)),
      list(c(-pi, pi, pi), c(-pi, pi, -pi)),
      list(c(-pi, 0, pi), c(-pi, 0, -pi)),
      list(c(pi, 0, pi), c(pi, 0, -pi)))
  else
    curves <- list()
  plot3D(radius=pi, points=xzxs, curves=curves, ...)
}

#' Breaking a curve into pieces that don't cross the boundary of a ball plot.
#'
#' Probably you don't want to use this. It's used internally in our plotting functions. Anyway, if two successive points are very far apart (more than 1.5 times the radius of the plot), then the curve is broken between them. Two antipodal boundary points are inserted, so that one piece of the curve ends on the boundary, and the next piece of the curve starts at the antipode. To clarify, if the original curve never needs breaking, then the returned list will have length 1.
#' @param vs A list of 3D real vectors.
#' @param radius A real number (positive). The radius of the plot.
#' @param ballFromRotation A function from the set of special orthogonal 3x3 real matrices to the set of 3D real vectors. A right inverse to rotationFromBall.
#' @param rotationFromBall A function from the set of 3D real vectors to the set of special orthogonal 3x3 real matrices. A left inverse to ballFromRotation.
#' @return A list of lists of 3D real vectors.
rotBallCurves <- function(vs, radius, ballFromRotation, rotationFromBall) {
  n <- length(vs)
  if (n <= 1)
    list(vs)
  else {
    curves <- list()
    curve <- list(vs[[1]])
    boundSquared <- (1.5 * radius)^2
    for (i in 1:(n - 1))
      if (sum((vs[[i]] - vs[[i + 1]])^2) <= boundSquared)
        # This segment of the curve does not need to be broken.
        curve[[length(curve) + 1]] <- vs[[i + 1]]
    else {
      # The curve must be broken between the ith and (i + 1)th elements.
      q <- rotationFromBall(vs[[i]])
      r <- rotationFromBall(vs[[i + 1]])
      u <- rotAxisAngleFromMatrix(r %*% t(q))[1:3]
      # Finding the boundary means finding t such that aa t^2 + bb t + cc == 0,
      # where aa, bb, cc have been derived using Mathematica.
      cc <- 1 + q[1, 1] + q[2, 2] + q[3, 3]
      bb <- u[1] * (q[2, 3] - q[3, 2]) + u[2] * (q[3, 1] - q[1, 3]) + u[3] * (q[1, 2] - q[2, 1])
      bb <- 2 * bb
      aa <- u[1] * u[2] * (q[1, 2] + q[2, 1]) - u[3]^2 * (q[1, 1] + q[2, 2])
      aa <- aa + u[1] * u[3] * (q[1, 3] + q[3, 1]) - u[2]^2 * (q[1, 1] + q[3, 3])
      aa <- aa + u[2] * u[3] * (q[2, 3] + q[3, 2]) - u[1]^2 * (q[2, 2] + q[3, 3])
      aa <- 2 * aa + cc
      sols <- realQuadraticSolutions(aa, bb, cc)
      if (length(sols) == 0 || length(sols) == 3) {
        # Degenerate case. Just break the segment here.
        curves[[length(curves) + 1]] <- curve
        curve <- list(vs[[i + 1]])
      } else {
        # Typical case. Find a ball vector at the boundary.
        s <- sols[[1]]
        a <- arcCos((1 - s^2) / (1 + s^2))
        r <- rotMatrixFromAxisAngle(c(u, a)) %*% q
        v <- ballFromRotation(r)
        # Determine how to pair v and -v with vs[[i]] and vs[[i + 1]].
        if (sum((vs[[i]] - v)^2) < sum((vs[[i]] + v)^2))
          w <- -v
        else {
          w <- v
          v <- -v
        }
        # Finish this curve and start the next.
        curve[[length(curve) + 1]] <- v
        curves[[length(curves) + 1]] <- curve
        curve <- list(w, vs[[i + 1]])
      }
    }
    curves[[length(curves) + 1]] <- curve
    curves
  }
}

#' Infrastucture for equal-volume and similar plots.
#'
#' Unless you are designing your own system of plotting, you probably do not want to call this function. Call rotEqualVolumePlot, etc. instead.
#' @param radius A real number (positive). If NULL, then a radius is chosen to contain all points, curves, and triangles.
#' @param points A list of 3D real vectors.
#' @param curves A list of lists of 3D real vectors.
#' @param triangles A list of length-3 lists of 3D real vectors.
#' @param colors A list of strings (colors). Used to color the points only. Colors are recycled, just as in rgl.points and rgl.spheres.
#' @param ballFromRotation A function from the set of rotation matrices to the set of 3D real vectors. A right inverse to rotationFromBall.
#' @param rotationFromBall A function from the set of 3D real vectors to the set of rotation matrices. A left inverse to ballFromRotation.
#' @param boundaryAlpha A real number (between 0 and 1). The opacity of the boundary sphere. A value of 0 turns off the sphere entirely.
#' @param simplePoints A logical. Whether to plot points as points or as spheres. Spheres are higher-quality and slower.
#' @param backgroundColor String (color). Color of background and fog.
#' @param curveColor A string (color). Color of curves.
#' @param curveWidth A real number (positive). Width of curves in pixels.
#' @param fogStyle A string, either "none", "linear", "exp", or "exp2". The style of fog used to suggest depth. See rgl.bg.
#' @param pointSize A real number (positive). The size of the points or spheres used to depict points. If simplePoints, then measured in pixels with default 3. If not simplePoints, then measured in the same units as the radius of the plot, with default 0.02 * radius.
#' @param trianglesRaw Either NULL or a vector of 9 * m real numbers. An alternative way to pass m triangles: All of the x-coordinates of vertices, then all of the ys, then all of the zs.
#' @param axesColors A character vector of length 3. The colors for the axes.
#' @return NULL.
rotBallPlot <- function(radius=NULL, points, curves, triangles, colors, ballFromRotation, rotationFromBall, boundaryAlpha=0.1, simplePoints=FALSE, backgroundColor="black", curveColor="white", curveWidth=1, fogStyle="linear", pointSize=NULL, trianglesRaw=NULL, axesColors=c("red", "green", "blue")) {
  # Initialize the window. The user will close it later.
  rgl.open()
  rgl.bg(color=backgroundColor, fogtype=fogStyle)
  # Figure out the radius and pointSize if necessary.
  if (is.null(radius)) {
    pointRadius <- 0
    curveRadius <- 0
    triangleRadius <- 0
    f <- function(pointList) {max(abs(simplify2array(pointList)))}
    if (length(points) >= 1)
      pointRadius <- f(points)
    if (length(curves) >= 1)
      curveRadius <- max(sapply(curves, f))
    if (length(triangles) >= 1)
      triangleRadius <- max(sapply(triangles, f))
    radius <- max(pointRadius, curveRadius, triangleRadius)
  }
  if (is.null(pointSize)) {
    if (simplePoints)
      pointSize <- 3
    else
      pointSize <- 0.02 * radius
  }
  # Draw the three coordinate axes.
  xs <- radius * c(0, 1, 1, 1, 1, 1, -1, -1, -1, -1)
  ys <- radius * c(0, 0, -0.1, 0.1, 0, 0, -0.1, 0.1, 0, 0)
  zs <- radius * c(0, 0, 0, 0, -0.1, 0.1, 0, 0, -0.1, 0.1)
  if (radius > 1) {
    xs <- c(xs, as.numeric(sapply(1:radius, function(x) c(x, x, x, x))))
    ys <- c(ys, replicate(floor(radius), c(-0.1, 0.1, 0, 0)))
    zs <- c(zs, replicate(floor(radius), c(0, 0, -0.1, 0.1)))
  }
  rgl.lines(x=xs, y=ys, z=zs, color=axesColors[[1]], lwd=3)
  rgl.lines(x=zs, y=xs, z=ys, color=axesColors[[2]], lwd=3)
  rgl.lines(x=ys, y=zs, z=xs, color=axesColors[[3]], lwd=3)
  # Draw the points.
  if (length(points) >= 1) {
    xs <- sapply(points, function(p) {p[1]})
    ys <- sapply(points, function(p) {p[2]})
    zs <- sapply(points, function(p) {p[3]})
    if (simplePoints)
      rgl.points(x=xs, y=ys, z=zs, color=colors, size=pointSize)
    else
      rgl.spheres(x=xs, y=ys, z=zs, radius=pointSize, color=colors, lit=FALSE)
  }
  # Draw the curves.
  if (length(curves) >= 1)
    for (curve in curves) {
      curvesNew <- rotBallCurves(curve, radius, ballFromRotation, rotationFromBall)
      for (cur in curvesNew) {
        xs <- sapply(cur, function(p) {p[1]})
        ys <- sapply(cur, function(p) {p[2]})
        zs <- sapply(cur, function(p) {p[3]})
        rgl.linestrips(x=xs, y=ys, z=zs, color=curveColor, lwd=curveWidth)
      }
    }
  # Draw the triangles.
  if (length(triangles) >= 1) {
    xs <- sapply(triangles, function(tri) {sapply(tri, function(p) {p[1]})})
    ys <- sapply(triangles, function(tri) {sapply(tri, function(p) {p[2]})})
    zs <- sapply(triangles, function(tri) {sapply(tri, function(p) {p[3]})})
    rgl.triangles(x=xs, y=ys, z=zs)#, alpha=1.0, back="cull")
  }
  if (!is.null(trianglesRaw)) {
    numVert <- length(trianglesRaw) / 3
    xs <- trianglesRaw[1:numVert]
    ys <- trianglesRaw[(numVert + 1):(2 * numVert)]
    zs <- trianglesRaw[(2 * numVert + 1):(3 * numVert)]
    rgl.triangles(x=xs, y=ys, z=zs)
  }
  # Draw the container sphere.
  if (boundaryAlpha > 0)
    rgl.spheres(x=0, y=0, z=0, radius=radius, alpha=boundaryAlpha, back="cull")
  NULL
}

#' Axis-angle plot of the space of rotations.
#'
#' The plot represents each rotation as a point, whose direction from the origin indicates the axis of rotation, and whose distance from the origin indicates the angle of rotation. Angles between 0 and pi are used. So SO(3) plots as a ball of radius pi. This plot does not have an equal-angle or equal-volume property. Curves that cross the boundary of the plot are clipped automatically. On the other hand, triangles are not clipped.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Plotting options: boundaryAlpha, simplePoints, etc. See rotBallPlot for details.
#' @return NULL.
rotAxisAnglePlot <- function(points=list(), curves=list(), triangles=list(), colors=c("white"), ...) {
  aaFromRotation <- function(r) rotLeftTangentFromMatrix(r, diag(c(1, 1, 1)))
  rotationFromAA <- function(v) rotMatrixFromLeftTangent(v, diag(c(1, 1, 1)))
  pointsNew <- lapply(points, aaFromRotation)
  curvesNew <- lapply(curves, function(curve) lapply(curve, aaFromRotation))
  trianglesNew <- lapply(triangles, function(tri) lapply(tri, aaFromRotation))
  rotBallPlot(pi, pointsNew, curvesNew, trianglesNew, colors, aaFromRotation, rotationFromAA, ...)
}

#' Equal-angle plot of the space of rotations.
#'
#' The plot represents each rotation as a point, whose direction from the origin indicates the axis of rotation, and whose distance from the origin indicates the tangent quarter-angle of rotation. All positive angles are potentially used. So SO(3) plots as a ball of radius 1.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Plotting options: boundaryAlpha, simplePoints, etc. See rotBallPlot for details.
#' @return NULL.
rotEqualAnglePlot <- function(points=list(), curves=list(), triangles=list(), colors=c("white"), ...) {
  pointsNew <- lapply(points, rotEqualAngleFromMatrix)
  curvesNew <- lapply(curves, function(curve) lapply(curve, rotEqualAngleFromMatrix))
  trianglesNew <- lapply(triangles, function(tri) lapply(tri, rotEqualAngleFromMatrix))
  rotBallPlot(1, pointsNew, curvesNew, trianglesNew, colors, rotEqualAngleFromMatrix, rotMatrixFromEqualAngle, ...)
}

#' Rodrigues plot of the space of rotations.
#'
#' The plot represents each rotation as a point, whose direction from the origin indicates the axis of rotation, and whose distance from the origin indicates the tangent half-angle of rotation. All positive angles are potentially used. So SO(3) plots as a ball of infinite radius.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Plotting options: boundaryAlpha, simplePoints, etc. See rotBallPlot for details.
#' @return NULL.
rotRodriguesPlot <- function(points=list(), curves=list(), triangles=list(), colors=c("white"), ...) {
  pointsNew <- lapply(points, rotRodriguesFromMatrix)
  curvesNew <- lapply(curves, function(curve) lapply(curve, rotRodriguesFromMatrix))
  trianglesNew <- lapply(triangles, function(tri) lapply(tri, rotRodriguesFromMatrix))
  rotBallPlot(NULL, pointsNew, curvesNew, trianglesNew, colors, rotRodriguesFromMatrix, rotMatrixFromRodrigues, ...)
}

#' Equal-volume plot of the space of rotations, presented in equal-volume coordinates.
#'
#' Probably you don't want to use this. Look at rotEqualVolumePlot first. Each rotation is passed to this function as an equal-volume vector, as produced by rotEqualVolumeFromMatrix, for example. The function simply plots those vectors, in a ball of radius approximately 0.62. The plot is equal-volume; it accurately represents volumes of regions of SO(3). Curves that cross the boundary of the plot are clipped automatically. On the other hand, triangles are not clipped.
#' @param points A list of 3D real vectors.
#' @param curves A list of lists of 3D real vectors.
#' @param triangles A list of length-3 lists of 3D real vectors.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Plotting options: boundaryAlpha, simplePoints, etc. See rotBallPlot for details.
#' @return NULL.
rotNativeEqualVolumePlot <- function(points=list(), curves=list(), triangles=list(), colors=c("white"), ...) {
  rotBallPlot(rotEqualVolumeRadius, points, curves, triangles, colors, rotEqualVolumeFromMatrix, rotMatrixFromEqualVolume, ...)
}

#' Equal-volume plot of the space of rotations.
#'
#' Each rotation is passed to this function as a special orthogonal matrix. The function converts to volumetric representation and then simply calls volumetricPlot.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Plotting options: boundaryAlpha, simplePoints, etc. See rotBallPlot for details.
#' @return NULL.
rotEqualVolumePlot <- function(points=list(), curves=list(), triangles=list(), colors=c("white"), ...) {
  pointsNew <- lapply(points, rotEqualVolumeFromMatrix)
  curvesNew <- lapply(curves, function(curve) {lapply(curve, rotEqualVolumeFromMatrix)})
  trianglesNew <- lapply(triangles, function(tri) {lapply(tri, rotEqualVolumeFromMatrix)})
  rotNativeEqualVolumePlot(points=pointsNew, curves=curvesNew, triangles=trianglesNew, colors, ...)
}



### HIGHER-LEVEL PLOTS ###

#' Convenience shortcut for plotting two sets of rotations.
#' 
#' @param pointsA A list of rotations.
#' @param pointsB A list of rotations.
#' @param colorA Character. A color, of the sort used in all R graphics routines.
#' @param colorB Character. A color, of the sort used in all R graphics routines.
#' @param ... Other arguments to be passed to the underlying rotEqualVolumePlot.
#' @return NULL.
rotEqualVolumePlotTwo <- function(pointsA, pointsB, colorA="white", colorB="white", ...) {
  rotEqualVolumePlot(c(pointsA, pointsB), 
                     colors=c(replicate(length(pointsA), colorA), replicate(length(pointsB), colorB)),
                     ...)
}

#' Convenience shortcut for plotting three sets of rotations.
#' 
#' @param pointsA A list of rotations.
#' @param pointsB A list of rotations.
#' @param colorA Character. A color, of the sort used in all R graphics routines.
#' @param colorB Character. A color, of the sort used in all R graphics routines.
#' @param colorC Character. A color, of the sort used in all R graphics routines.
#' @param ... Other arguments to be passed to the underlying rotEqualVolumePlot.
#' @return NULL.
rotEqualVolumePlotThree <- function(pointsA, pointsB, pointsC, colorA="white", colorB="white", colorC="white", ...) {
  rotEqualVolumePlot(c(pointsA, pointsB, pointsC),
                     colors=c(replicate(length(pointsA), colorA), 
                              replicate(length(pointsB), colorB), 
                              replicate(length(pointsB), colorC)),
                     ...)
}

#' Equal-angle plot of geodesic regression results.
#' 
#' The user may choose to draw extra spurs, joining the Rs to their corresponding predictions on the curve.
#' @param xs A vector of real numbers.
#' @param rs A list of rotation matrices.
#' @param regr The output of rotGeodesicRegression with those xs and rs.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param spurs A list of rotation matrices. When not NULL, it's usually equal to Rs.
#' @param ... Other plotting options to be passed to rotEqualAnglePlot.
#' @return NULL.
rotGeodesicRegressionPlot <- function(xs, rs, regr, colors=hues(xs), spurs=NULL, ...) {
  rotA <- rotExp(min(xs) * regr$m) %*% regr$b
  rotB <- rotExp(mean(range(xs)) * regr$m) %*% regr$b
  rotC <- rotExp(max(xs) * regr$m) %*% regr$b
  curves <- list(rotGeodesicPoints(rotA, rotB, 20), rotGeodesicPoints(rotB, rotC, 20))
  if (!is.null(spurs)) {
    spurCurves <- lapply(1:length(xs), function(i) {rotGeodesicPoints(spurs[[i]], regr$prediction(xs[[i]]), 10)})
    curves <- c(curves, spurCurves)
  }
  rotEqualAnglePlot(points=rs, curves=curves, colors=colors, ...)
}

#' Equal-angle plot of Mahalanobis inference results.
#' 
#' Warning: Ellipsoid does not clip to the boundary of the plot correctly.
#' @param points A list of rotation matrices.
#' @param center A rotation matrix. Typically the $center from the inference.
#' @param leftCovarInv A 3x3 real matrix (symmetric, positive-definite). Typically the $leftCovarInv from the inference.
#' @param level A real number. Typically $q095^2 from the inference.
#' @param numNonAdapt A real number (non-negative integer). The number of refinements to the sphere that is deformed into the ellipsoid. Incrementing numNonAdapt improves visual quality but increases time and memory requirements by a factor of four.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Additional plotting options to be passed to rotEqualAnglePlot.
#' @return NULL.
rotEllipsoidPlot <- function(points, center, leftCovarInv, level, numNonAdapt=3, colors=c("white"), ...) {
  tris <- rotEllipsoidTriangles(center, leftCovarInv, level, numNonAdapt=numNonAdapt)
  rotEqualAnglePlot(points=points, triangles=tris, colors=colors, ...)
}

# Helper function for rotEqualAreaTickPlot.
rotEqualAreaTick <- function(v, tickSize) {
  r <- sqrt(v[[1]]^2 + v[[2]]^2)
  if (r <= 0 || 1 - r^2 <= 0)
    tick <- c(0, 0, 0)
  else {
    scalar <- sqrt(1 - r^2) / r
    tick <- c(scalar * v[[1]], scalar * v[[2]], v[[3]] / -scalar)
  }
  list(lower(v), rayNormalized(lower(v) + tickSize * tick))
}

#' Equal-area plot, attempting to show all degrees of freedom.
#' 
#' In each rotation, regards the first row as +-down-pole and the third row as +-ray. Depicts each rotation as a great circle (perpendicular to the pole), a point on that great circle (corresponding to the ray), and a tick mark (indicating the direction of the ray, either toward the center if up or away if down). Really intended for representing faults-with-slip in geology, where the second row is then the vorticity of fault slip. But could be useful in more abstract rotations?
#' @param rs A list of rotation matrices.
#' @param numSteps A real number (positive integer). The number of steps used in drawing each great circle.
#' @param tickSize A real number (non-negative). The length of the ticks.
#' @param ... Additional options to be passed to the underlying lineEqualAreaPlot.
#' @return NULL.
rotEqualAreaTickPlot <- function(rs, numSteps=60, tickSize=0.1, ...) {
  greats <- lapply(rs, function(r) rayGreatCircle(r[1,], numSteps))
  hangs <- lapply(rs, function(r) {if (r[[1, 3]] < 0) r[3,] else -r[3,]})
  ticks <- lapply(hangs, rotEqualAreaTick, tickSize)
  lineEqualAreaPlot(points=hangs, curves=c(greats, ticks), ...)
}



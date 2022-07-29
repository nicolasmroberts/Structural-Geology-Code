


# Copyright 2016-2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### ELLIPSOID PLOTS ###

#' Visualization of ellipsoids, as surfaces with curves drawn on them.
#' 
#' @param rots A list of 3x3 real matrices (special orthogonal). The ellipsoid orientations, as given by the $rotation field of an ellipsoid.
#' @param as A list of 3D real vectors. The ellipsoid semi-axis lengths, in order corresponding to the rows of the rots, as given by the $a field of an ellipsoid.
#' @param centers A list of 3D real vectors. The locations of the ellipsoid centers.
#' @param numNonAdapt A real number (non-negative integer). The number of refinements to use. Each refinement makes the ellipsoids smoother, but increases time and memory requirements by a factor of four.
#' @param ... Other arguments to be passed to the underlying plot3D.
#' @return NULL.
ellEllipsoidPlot <- function(rots, as, centers=replicate(length(rots), c(0, 0, 0), simplify=FALSE), numNonAdapt=4, ...) {
  # This function applies the ith ellipsoid's transformation to the given vector v.
  transf <- function(v, i) {centers[[i]] + as.numeric(t(rots[[i]]) %*% (as[[i]] * v))}
  # Make surfaces.
  sphere <- rayTetrahedralSphere(numNonAdapt)
  triangles <- unlist(
    lapply(1:length(as), function(i) lapply(sphere, function(tri) lapply(tri, transf, i))),
    recursive=FALSE, use.names=FALSE)
  # Make curves.
  thetas <- seq(from=0, to=(2 * pi), length.out=100)
  circle12 <- lapply(thetas, function(theta) c(cos(theta), sin(theta), 0))
  circle23 <- lapply(thetas, function(theta) c(0, cos(theta), sin(theta)))
  circle31 <- lapply(thetas, function(theta) c(sin(theta), 0, cos(theta)))
  circles <- list(circle12, circle23, circle31)
  curves <- unlist(
    lapply(1:length(as), function(i) lapply(circles, function(circle) lapply(circle, transf, i))),
    recursive=FALSE, use.names=FALSE)
  plot3D(triangles=triangles, curves=curves, ...)
}

#' Equal-area plot of ellipsoid axes.
#' 
#' Short axes are shown as circles, intermediate as triangles, long as squares. Warning: Curves are not well tested.
#' @param rots A list of 3x3 real matrices (special orthogonal). The ellipsoid orientations, as given by the $rotation field of an ellipsoid.
#' @param as A list of 3D real vectors. The ellipsoid semi-axis lengths, in order corresponding to the rows of the rots, as given by the $a field of an ellipsoid. Alternatively, you can pass logA; this information is used only to determine the order of the axes.
#' @param rotCurves A list of lists of 3x3 real matrices (special orthogonal). Like rots, but curves rather than points.
#' @param aCurves A list of lists of 3D real vectors. Like as, but curves rather than points.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @return NULL.
ellEqualAreaPlot <- function(rots, as, rotCurves=list(), aCurves=list(), colors=c("black")) {
  # Prepare to plot points based on rots and as.
  f <- function(i, rs, as) {
    ord <- order(as[[i]])
    list(rs[[i]][ord[[1]],], rs[[i]][ord[[2]],], rs[[i]][ord[[3]],])
  }
  points <- unlist(lapply(1:length(rots), f, rots, as), recursive=FALSE, use.names=FALSE)
  # Prepare to plot curves based on curvesRots and curvesAs.
  if (length(rotCurves) >= 1 && length(rotCurves) == length(aCurves)) {
    curves <- lapply(1:length(rotCurves), function(j) lapply(1:length(rotCurves[[j]]), f, rotCurves[[j]], aCurves[[j]]))
    curves1 <- lapply(curves, function(curve) lapply(curve, function(tri) tri[[1]]))
    curves2 <- lapply(curves, function(curve) lapply(curve, function(tri) tri[[2]]))
    curves3 <- lapply(curves, function(curve) lapply(curve, function(tri) tri[[3]]))
    curves <- c(curves1, curves2, curves3)
  } else
    curves <- list()
  # Plot.
  newColors <- as.character(sapply(colors, function(s) c(s, s, s)))
  lineEqualAreaPlot(points, curves=curves, colors=newColors, shapes=c("c", "t", "s"))
}

#' Equal-volume plot of ellipsoid orientations.
#' 
#' @param rots A list of 3x3 real matrices (special orthogonal). The ellipsoid orientations, as given by the $rotation field of an ellipsoid.
#' @param as A list of 3D real vectors. The ellipsoid semi-axis lengths, in order corresponding to the rows of the rots, as given by the $a field of an ellipsoid. Alternatively, you can pass logA; this information is used only to determine the order of the axes.
#' @param rotCurves A list of lists of 3x3 real matrices (special orthogonal). Like rots, but curves rather than points.
#' @param aCurves A list of lists of 3D real vectors. Like as, but curves rather than points.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @param ... Other arguments to be passed to the underlying oriEqualVolumePlot.
#' @return NULL.
ellEqualVolumePlot <- function(rots, as, rotCurves=list(), aCurves=list(), colors=c("white"), ...) {
  # The rotations are permuted into this row-order: short, long, intermediate. To match other plane-line stuff in our library.
  f <- function(i, rs, as) {
    ord <- order(as[[i]])
    short <- rs[[i]][ord[[1]],]
    long <- rs[[i]][ord[[3]],]
    rbind(short, long, cross(short, long))
  }
  points <- lapply(1:length(rots), f, rots, as)
  if (length(rotCurves) >= 1 && length(rotCurves) == length(aCurves)) {
    curves <- lapply(1:length(rotCurves), function(j) lapply(1:length(rotCurves[[j]]), f, rotCurves[[j]], aCurves[[j]]))
    f <- function(curve) {
      cur <- list(curve[[1]])
      for (r in curve[2:length(curve)])
        cur[[length(cur) + 1]] <- oriNearestRepresentative(r, cur[[length(cur)]], oriLineInPlaneGroup)
      cur
    }
    curves <- lapply(curves, f)
  } else
    curves <- list()
  oriEqualVolumePlot(points=points, curves=curves, group=oriLineInPlaneGroup, colors=colors, ...)
}

#' A distorted version of the Hsu-Nadai plot of ellipsoid shapes.
#' 
#' This is a polar plot, in which the radial coordinate is octahedral shear strain and the angular coordinate is Lode's parameter. This plot is similar to, but not identical to, the Hsu-Nadai plot. See ellHsuNadaiPlot for the real thing.
#' @param logAs A list of 3D real vectors. The ellipsoid semi-axis log-lengths.
#' @param curves A list of lists of 3D real vectors. Like logAs, but curves rather than points.
#' @param es A real number (positive) or NULL. If a number, then that is the radius of the plot in the E_s direction. If NULL, then the radius of the plot is inferred from the points (not the curves, currently).
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @return NULL.
ellWrongHsuNadaiPlot <- function(logAs, curves=list(), es=NULL, colors=c("black")) {
  ess <- sapply(logAs, ellOctahedralShearStrain)
  nus <- sapply(logAs, ellLodeNu)
  esCurves <- lapply(curves, function(curve) sapply(curve, ellOctahedralShearStrain))
  nuCurves <- lapply(curves, function(curve) sapply(curve, ellLodeNu))
  # Make the plot window.
  if (is.null(es))
    es <- max(c(ess, 1))
  plot.new()
  plot.window(xlim=c(-0.55 * es, 0.55 * es), ylim=c(-0.05 * es, 1.05 * es))
  # Plot the points.
  if (length(logAs) >= 1) {
    xs <- sapply(1:length(logAs), function(i) ess[[i]] * cos(pi / 2 - nus[[i]] * pi / 6))
    ys <- sapply(1:length(logAs), function(i) ess[[i]] * sin(pi / 2 - nus[[i]] * pi / 6))
    points(xs, ys, col=colors, pch=c(19))
  }
  # Plot the curves.
  if (length(curves) >= 1)
    for (j in 1:length(curves)) {
      xs <- sapply(1:length(curves[[j]]), function(i) esCurves[[j]][[i]] * cos(pi / 2 - nuCurves[[j]][[i]] * pi / 6))
      ys <- sapply(1:length(curves[[j]]), function(i) esCurves[[j]][[i]] * sin(pi / 2 - nuCurves[[j]][[i]] * pi / 6))
      lines(xs, ys)
    }
  # Plot the boundary.
  xys <- sapply(0:30, function(t) {
    theta <- (t / 30) * (pi / 3) + (pi / 3)
    es * c(cos(theta), sin(theta))
  })
  lines(xys[1,], xys[2,])
  lines(c(-0.5 * es, 0, 0.5 * es), c(sqrt(3) / 2 * es, 0, sqrt(3) / 2 * es))
  # Plot some tick marks.
  if (es >= 1)
    for (i in 1:es) {
      xs <- c(i * 0.5, i * 0.5 + sqrt(3) * 0.5 * 0.05)
      ys <- c(i * sqrt(3) * 0.5, i * sqrt(3) * 0.5 - 0.5 * 0.05)
      lines(xs, ys)
      lines(-xs, ys)
    }
}

#' Hsu-Nadai plot of ellipsoid shapes.
#' 
#' @param logAs A list of 3D real vectors. The ellipsoid semi-axis log-lengths.
#' @param curves A list of lists of 3D real vectors. Like logAs, but curves rather than points.
#' @param es A real number (positive) or NULL. If a number, then that is the radius of the plot in the E_s direction. If NULL, then the radius of the plot is inferred from the points (not the curves, currently).
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @param shapes A pch-style argument.
#' @return NULL.
ellHsuNadaiPlot <- function(logAs, curves=list(), es=NULL, colors=c("black"), shapes=c(19)) {
  x <- function(logA) {-sqrt(1.5) * max(logA - sum(logA) / 3) - sqrt(1.5) * min(logA - sum(logA) / 3)}
  y <- function(logA) {sqrt(0.5) * max(logA - sum(logA) / 3) - sqrt(0.5) * min(logA - sum(logA) / 3)}
  # Make the plot window.
  if (is.null(es))
    es <- max(c(1, sapply(logAs, ellOctahedralShearStrain)))
  plot.new()
  plot.window(xlim=c(-0.55 * es, 0.55 * es), ylim=c(-0.05 * es, 1.05 * es))
  # Plot the points.
  if (length(logAs) >= 1)
    points(sapply(logAs, x), sapply(logAs, y), col=colors, pch=shapes)
  # Plot the curves.
  if (length(curves) >= 1)
    for (j in 1:length(curves))
      lines(sapply(curves[[j]], x), sapply(curves[[j]], y))
  # Plot the boundary.
  xys <- sapply(0:30, function(t) {
    theta <- (t / 30) * (pi / 3) + (pi / 3)
    es * c(cos(theta), sin(theta))
  })
  lines(xys[1,], xys[2,])
  lines(c(-0.5 * es, 0, 0.5 * es), c(sqrt(3) / 2 * es, 0, sqrt(3) / 2 * es))
  # Plot some tick marks.
  if (es >= 1)
    for (i in 1:es) {
      xs <- c(i * 0.5, i * 0.5 + sqrt(3) * 0.5 * 0.05)
      ys <- c(i * sqrt(3) * 0.5, i * sqrt(3) * 0.5 - 0.5 * 0.05)
      lines(xs, ys)
      lines(-xs, ys)
    }
}

#' Hsu-Nadai plot of ellipsoid shapes, with a third dimension specified by the user.
#' 
#' @param logAs A list of 3D real vectors. The ellipsoid semi-axis log-lengths.
#' @param zs A vector of real numbers. The coordinates of the points in the direction perpendicular to the Hsu-Nadai plot.
#' @param es A real number (positive) or NULL. If a number, then that is the radius of the plot in the E_s direction. If NULL, then the radius of the plot is inferred from the points.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @param ... Other arguments to pass to the underlying plot3D.
#' @return NULL.
ellHsuNadaiScalarPlot <- function(logAs, zs, es=NULL, colors=c("white"), ...) {
  x <- function(logA) {-sqrt(1.5) * max(logA - sum(logA) / 3) - sqrt(1.5) * min(logA - sum(logA) / 3)}
  y <- function(logA) {sqrt(0.5) * max(logA - sum(logA) / 3) - sqrt(0.5) * min(logA - sum(logA) / 3)}
  # Determine how big the plot will be.
  if (is.null(es))
    es <- max(c(1, sapply(logAs, ellOctahedralShearStrain)))
  z <- max(abs(zs))
  radius <- max(es, z)
  # Build the points.
  points <- lapply(1:length(logAs), function(i) c(x(logAs[[i]]), y(logAs[[i]]), zs[[i]]))
  # Build the curves.
  f <- function(t, zz) {
    theta <- (t / 30) * (pi / 3) + (pi / 3)
    c(es * c(cos(theta), sin(theta)), zz)
  }
  bottom <- lapply(0:30, f, -z)
  bottom <- c(bottom, list(c(-0.5 * es, sqrt(3) / 2 * es, -z), c(0, 0, -z), c(0.5 * es, sqrt(3) / 2 * es, -z)))
  middle <- lapply(0:30, f, 0)
  middle <- c(middle, list(c(-0.5 * es, sqrt(3) / 2 * es, 0), c(0, 0, 0), c(0.5 * es, sqrt(3) / 2 * es, 0)))
  top <- lapply(0:30, f, z)
  top <- c(top, list(c(-0.5 * es, sqrt(3) / 2 * es, z), c(0, 0, z), c(0.5 * es, sqrt(3) / 2 * es, z)))
  plot3D(radius=radius, points=points, curves=list(bottom, middle, top), colors=colors, ...)
}

#' Flinn plot of ellipsoid shapes.
#' 
#' @param as A list of 3D real vectors, with all entries positive. The ellipsoid semi-axis lengths.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @return NULL.
ellFlinnPlot <- function(as, colors=c("black")) {
  xs <- sapply(as, function(a) {(sum(a) - min(a) - max(a)) / min(a)})
  ys <- sapply(as, function(a) {max(a) / (sum(a) - min(a) - max(a))})
  plot(x=xs, y=ys, xlim=c(1, max(xs)), ylim=c(1, max(ys)),
       xlab="intermediate / short", ylab="long / intermediate")
}

#' Logarithmic Flinn plot (Ramsay plot) of ellipsoid shapes.
#' 
#' @param logAs A list of 3D real vectors. The ellipsoid semi-axis log-lengths.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @return NULL.
ellLogFlinnPlot <- function(logAs, colors=c("black")) {
  xs <- sapply(logAs, function(logA) {-max(logA - sum(logA) / 3) - 2 * min(logA - sum(logA) / 3)})
  ys <- sapply(logAs, function(logA) {2 * max(logA - sum(logA) / 3) + min(logA - sum(logA) / 3)})
  plot(x=xs, y=ys, xlim=c(0, max(xs)), ylim=c(0, max(ys)),
       xlab="log(intermediate / short)", ylab="log(long / intermediate)")
}

#' Jelinek plot of ellipsoid shapes.
#' 
#' This is a rectangular plot of Jelinek's Pj vs. Lode's nu. Tauxe (2010) called it the Jelinek plot, after Jelinek (1981).
#' @param logAs A list of 3D real vectors. The ellipsoid semi-axis log-lengths.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @return NULL.
ellJelinekPlot <- function(logAs, colors=c("black")) {
  plot(x=sapply(logAs, ellJelinekP), 
       y=sapply(logAs, ellLodeNu), 
       col=colors, xlab="P_j", ylab="nu", ylim=c(-1, 1))
}

#' Pair plot of ellipsoid vectors.
#' 
#' @param points A list of 5D or 6D real vectors. Ellipsoid vectors, as in the $vector field of an ellipsoid.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @param ... Other parameters to be passed to the underlying pairs function.
#' @param NULL.
ellPairsPlot <- function(points, colors=c("black"), ...) {
  pairs(t(simplify2array(points)), labels=c("v_1", "v_2", "v_3", "v_4", "v_5"), col=colors, ...)
}

#' Convenience function for quick plot of pairs of principal component weights.
#' 
#' @param points A list of 5D or 6D real vectors. Ellipsoid vectors, as in the $vector field of an ellipsoid.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @param ... Other parameters to be passed to the underlying pairs function.
#' @param NULL.
ellPCAPairsPlot <- function(points, colors=c("black"), ...) {
  pca <- ellPrincipalComponentAnalysis(points)
  pcs <- lapply(points, function(v) as.numeric(t(pca$rotation) %*% (v - pca$center)))
  pairs(t(simplify2array(pcs)), labels=c("pc_1", "pc_2", "pc_3", "pc_4", "pc_5"), col=colors, ...)
}

#' 2D or 3D plot of ellipsoid vectors.
#' 
#' The 2D case is like a single panel of ellPairsPlot. Warning: In 2D, curves are not well tested.
#' @param ijk A 2D or 3D vector of real numbers (positive integers). These should be in 1, ..., d, where d is the dimension of the vectors. They select out which coordinates of the vectors to display. For example, c(1, 2, 3) indicates to make a 3D plot of the first three vector coordinates.
#' @param points A list of 5D or 6D real vectors. Ellipsoid vectors, as in the $vector field of an ellipsoid.
#' @param curves A list of lists of 5D or 6D real vectors. Like points, but curves.
#' @param colors Character or NULL. A vector of colors for coloring the points, as in all R graphics functions. If NULL, then defaults to black in 2D or white in 3D.
#' @param ... Other parameters to be passed to the underlying plot or plot3D function.
#' @param NULL.
ellVectorPlot <- function(ijk, points=list(), curves=list(), colors=NULL, ...) {
  pointsNew <- lapply(points, function(v) v[ijk])
  curvesNew <- lapply(curves, function(curve) lapply(curve, function(v) v[ijk]))
  if (length(ijk) == 3) {
    if (is.null(colors))
      colors="white"
    plot3D(points=pointsNew, curves=curvesNew, colors=colors, ...)
  } else {
    if (is.null(colors))
      colors="black"
    plot(t(simplify2array(pointsNew)), col=colors,
         xlab=paste0("ellipsoid v_", as.character(ijk[[1]])),
         ylab=paste0("ellipsoid v_", as.character(ijk[[2]])), ...)
    for (curve in curvesNew)
      lines(t(simplify2array(curve)))
  }
}

#' 2D or 3D plot of ellipsoid vectors, as weights after principal component analysis.
#' 
#' The 2D case is like a single panel of ellPCAPairsPlot. Warning: In 2D, curves are not well tested.
#' @param ijk A 2D or 3D vector of real numbers (positive integers). These should be in 1, ..., d, where d is the dimension of the vectors. They select out which coordinates of the vectors to display. For example, c(1, 2, 3) indicates to make a 3D plot of the first three vector coordinates.
#' @param points A list of 5D or 6D real vectors. Ellipsoid vectors, as in the $vector field of an ellipsoid.
#' @param curves A list of lists of 5D or 6D real vectors. Like points, but curves.
#' @param colors Character or NULL. A vector of colors for coloring the points, as in all R graphics functions. If NULL, then defaults to black in 2D or white in 3D.
#' @param ... Other parameters to be passed to the underlying plot3D function. Ignored in the 2D case.
#' @param NULL.
ellPCAPlot <- function(ijk, points=list(), curves=list(), colors=NULL, ...) {
  pca <- ellPrincipalComponentAnalysis(points)
  transf <- function(v) as.numeric(t(pca$rotation) %*% (v - pca$center))
  pointsNew <- lapply(points, transf)
  curvesNew <- lapply(curves, function(curve) lapply(curve, transf))
  ellVectorPlot(ijk, points=pointsNew, curves=curvesNew, colors=colors, ...)
}






# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### INFRASTRUCTURE FOR PLOTTING FUNCTIONS ###

# Given 0 <= s <= 1, returns the line that is fraction s of the way from u to v. If u and v are perpendicular, then the ambiguity is resolved by treating them like rays.
lineInterpolation <- function(u, v, s=0.5) {
  if (dot(u, v) < 0)
    rayInterpolation(u, -v, s)
  else
    rayInterpolation(u, v, s)
}

# Returns a list of five boxes. Each box is a list of four vertices, in counter-clockwise order when viewed from above the sphere. Each vertex is a vector of four numbers: a unit 3D vector, with the value of f at that vector appended.
lineFunctionPlotBase <- function(f, ...) {
  a <- sqrt(2) / 2
  us <- list(c(a, 0, -a), c(0, a, -a), c(-a, 0, -a), c(0, -a, -a), c(1, 0, 0), c(0, 1, 0), c(-1, 0, 0), c(0, -1, 0))
  us <- lapply(us, function(u) c(u, f(u, ...)))
  list(
    list(us[[1]], us[[2]], us[[3]], us[[4]]),
    list(us[[1]], us[[5]], us[[6]], us[[2]]),
    list(us[[2]], us[[6]], us[[7]], us[[3]]),
    list(us[[3]], us[[7]], us[[8]], us[[4]]),
    list(us[[4]], us[[8]], us[[5]], us[[1]]))
}

# Given one box, returns a list of four boxes.
lineFunctionPlotRefinement <- function(box, f, ...) {
  mids <- lapply(1:4, function(i) lineInterpolation(box[[i]][1:3], box[[i %% 4 + 1]][1:3]))
  center <- lineInterpolation(lineInterpolation(mids[[1]], mids[[3]]), lineInterpolation(mids[[2]], mids[[4]]))
  mids <- lapply(mids, function(u) c(u, f(u, ...)))
  center <- c(center, f(center, ...))
  list(
    list(box[[1]], mids[[1]], center, mids[[4]]),
    list(box[[2]], mids[[2]], center, mids[[1]]),
    list(box[[3]], mids[[3]], center, mids[[2]]),
    list(box[[4]], mids[[4]], center, mids[[3]]))
}

# Given one box, returns a list of lines (pairs of unit vectors).
lineFunctionPlotCuts <- function(box, level) {
  cuts <- list()
  for (i in 1:4) {
    a <- box[[i]]
    b <- box[[i %% 4 + 1]]
    if (a[[4]] == level)
      cuts <- c(cuts, list(a[1:3]))
    else if (a[[4]] < level && b[[4]] > level)
      cuts <- c(cuts, list(lineInterpolation(a[1:3], b[1:3], (level - a[[4]]) / (b[[4]] - a[[4]]))))
    else if (b[[4]] < level && a[[4]] > level)
      cuts <- c(cuts, list(lineInterpolation(b[1:3], a[1:3], (level - b[[4]]) / (a[[4]] - b[[4]]))))
  }
  allPairs(cuts)
}

# Returns a list of length equal to that of levels. Each item in the list is a list of lines. Each line is a pair of unit vectors. Occasionally little 'stars' (complete graphs) may appear along a contour. If so, the non-adaptive refinement was not fine enough to resolve the contour in that box of the plot. So increase numNonAdapt.
lineFunctionPlotLineSets <- function(f, levels, numNonAdapt, ...) {
  # Non-adaptive refinement.
  boxes <- lineFunctionPlotBase(f, ...)
  if (numNonAdapt >= 1)
    for (i in 1:numNonAdapt) 
      boxes <- unlist(lapply(boxes, lineFunctionPlotRefinement, f, ...), recursive=FALSE)
    # Generate lines for this level of f.
    lineSets <- lapply(levels, function(level) unlist(lapply(boxes, lineFunctionPlotCuts, level), recursive=FALSE))
    lineSets
}

# lineFunctionPlotLineSets customized to the problem of Kamb contouring.
lineKambLineSets <- function(points, multiples=c(3, 6, 9, 12), k=3, degree=3, numNonAdapt=4) {
  # Compute the basic Kamb parameters.
  n <- length(points)
  p <- k^2 / (n + k^2)
  sigma <- squareRoot(n * p * (1 - p))
  r <- arcCos(1 - p)
  levels <- multiples * sigma
  # Compute the weighting function.
  if (degree == 0) {
    h <- function(p, u) {
      alpha <- lineDistance(p, u)
      if (alpha > r)
        0
      else
        1
    }
  } else {
    d0 <- 1 - cos(r)
    d1 <- -r * cos(r) + sin(r)
    d2 <- -2 + 2 * cos(r) - r^2 * cos(r) + 2 * r * sin(r)
    d3 <- 6 * r * cos(r) - r^3 * cos(r) - 6 * sin(r) + 3 * r^2 * sin(r)
    if (degree == 1)
      mat <- rbind(c(d0, d1, d2, d3), c(1, r, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1))
    else
      # Default to degree-3 weighting polynomial.
      mat <- rbind(c(d0, d1, d2, d3), c(0, 1, 0, 0), c(0, 1, 2 * r, 3 * r^2), c(1, r, r^2, r^3))
    coeffs <- solve(mat, c(1 - cos(r), 0, 0, 0))
    h <- function(p, u) {
      alpha <- lineDistance(p, u)
      if (alpha > r)
        0
      else
        coeffs[[1]] + coeffs[[2]] * alpha + coeffs[[3]] * alpha^2 + coeffs[[4]] * alpha^3
    }
  }
  # Contour-plot the function.
  f <- function(u) sum(sapply(points, h, u))
  lineFunctionPlotLineSets(f, levels, numNonAdapt)
}



### PLOTS THEMSELVES ###

#' Equal-area, lower-hemisphere plot of lines.
#' 
#' This function essentially sends everything to the lower hemisphere and then invokes rayEqualAreaPlot. See rayEqualAreaPlot for details.
#' @param points A list of rays.
#' @param curves A list of lists of rays.
#' @param colors A character vector.
#' @param shapes A character vector. Traditionally lineations are squares, intermediates are triangles, and poles to foliation are circles.
#' @return NULL
lineEqualAreaPlot <- function(points=list(), curves=list(), colors=c("black"), shapes=c("c")) {
  # Convert the chosen shapes into their underlying shape codes.
  f <- function(i) {
    j <- (i - 1) %% length(shapes) + 1
    if (shapes[[j]] == "s")
      15
    else if (shapes[[j]] == "t")
      17
    else if (shapes[[j]] == "c")
      19
    else
      46
  }
  underShapes <- sapply(1:length(points), f)
  # Send the points to the lower hemisphere.
  underPoints <- lapply(points, lower)
  # Break the curves and send them to the lower hemisphere.
  underCurves <- list()
  for (curve in curves) {
    curvesSigns <- rayCurvesUpperLower(curve)
    for (i in 1:length(curvesSigns$curves)) {
      if (curvesSigns$signs[[i]] == 1)
        underCurves[[length(underCurves) + 1]] <- lapply(curvesSigns$curves[[i]], function(v) -v)
      else
        underCurves[[length(underCurves) + 1]] <- curvesSigns$curves[[i]]
    }
  }
  plotEqualArea(points=underPoints, curves=underCurves, colors=colors, shapes=underShapes)
}

#' Equal-area, lower-hemisphere plot of lines, with a circle about each point.
#' 
#' This function essentially sends everything to the lower hemisphere and then invokes rayEqualAreaRadiusPlot.
#' @param points A list of lines, to be plotted as points.
#' @param radii A vector of real numbers. The radii of circles to be drawn about the points. Radii are measured in radians, along the surface of the unit sphere.
#' @param colors A character vector. See rayEqualAreaPlot.
#' @param shapes A character vector. See rayEqualAreaPlot.
#' @return NULL
lineEqualAreaRadiusPlot <- function(points=list(), radii=c(), colors=c("black"), shapes=c("c")) {
  lineEqualAreaPlot(points, colors=colors, curves=thread(raySmallCircle, points, radii))
}

#' Equal-area, lower-hemisphere plot of lines, with contours representing level sets of an arbitrary given function.
#' 
#' @param f An R function. Its input is a single line u, and optionally additional parameters '...'. Its output is a real number f(u, ...).
#' @param levels A vector of real numbers. The numbers c such that the contours f(u, ...) = c are plotted.
#' @param numNonAdapt A real number (non-negative integer). The number of non-adaptive refinements to make. These control the quality of the plot. Each increment to numNonAdapt makes the contours smoother, and helps resolve 'star defects' that sometimes appear along the contours. On the other hand, each increment increases the time and memory requirements by a factor of four.
#' @param points A list of lines, to be plotted as points. Need not have anything to do with f.
#' @param colors A character vector. See rayEqualAreaPlot.
#' @param shapes A character vector. See rayEqualAreaPlot.
#' @param ... Additional arguments to be passed to f, other than u.
#' @return NULL
lineEqualAreaFunctionPlot <- function(f, levels, numNonAdapt=4, points=list(), colors=c("black"), shapes=c("c"), ...) {
  lineSets <- lineFunctionPlotLineSets(f, levels, numNonAdapt, ...)
  lineEqualAreaPlot(points=points, curves=unlist(lineSets, recursive=FALSE), colors=colors, shapes=shapes)
}

#' Equal-area, lower-hemisphere plot of lines, with Kamb contours representing density.
#' 
#' @param points A list of lines, to be used in the computation of the Kamb contours.
#' @param multiples A vector of real numbers (positive). The multiples of sigma to be plotted as contours.
#' @param k A real number (positive). The arbitrary smoothing factor of Kamb (1959). I see no practical reason to mess with this, other than curiosity.
#' @param degree A real number (0, 1, or 3). The degree of the weighting polynomial. Higher degrees tend to produce smoother plots. I see no practical reason to mess with this, other than curiosity.
#' @param numNonAdapt A real number (non-negative integer). The number of non-adaptive refinements to make. These control the quality of the plot. Each increment to numNonAdapt makes the contours smoother, and helps resolve 'star defects' that sometimes appear along the contours. On the other hand, each increment increases the time and memory requirements by a factor of four.
#' @param colors A character vector. See rayEqualAreaPlot.
#' @param shapes A character vector. See rayEqualAreaPlot.
#' @param pointsToShow A list of lines, to be plotted as points.
#' @return NULL
lineKambPlot <- function(points, multiples=c(3, 6, 9, 12), k=3, degree=3, numNonAdapt=4, colors=c("black"), shapes=c("c"), pointsToShow=points) {
  lineSets <- lineKambLineSets(points=points, multiples=multiples, k=k, degree=degree, numNonAdapt=numNonAdapt)
  lineEqualAreaPlot(points=pointsToShow, curves=unlist(lineSets, recursive=FALSE), colors=colors, shapes=shapes)
}

#' Convenience shortcut for plotting two sets of lines.
#' 
#' @param pointsA A list of lines.
#' @param pointsB A list of lines.
#' @param colorA Character. A color, of the sort used in all R graphics routines.
#' @param colorB Character. A color, of the sort used in all R graphics routines.
#' @param shapeA Character. Either "c" for circle, "t" for triangle, or "s" for square.
#' @param shapeB Character. See shapeA.
#' @param curves A list of lists of lines. Curves to be plotted.
#' @return NULL.
lineEqualAreaPlotTwo <- function(pointsA, pointsB, colorA="black", colorB="black", shapeA="c", shapeB="c", curves=list()) {
  lineEqualAreaPlot(points=c(pointsA, pointsB), curves=curves,
                    shapes=c(replicate(length(pointsA), shapeA), replicate(length(pointsB), shapeB)),
                    colors=c(replicate(length(pointsA), colorA), replicate(length(pointsB), colorB)))
}

#' Convenience shortcut for plotting three sets of lines.
#' 
#' @param pointsA A list of lines.
#' @param pointsB A list of lines.
#' @param pointsC A list of lines.
#' @param colorA Character. A color, of the sort used in all R graphics routines.
#' @param colorB Character. A color, of the sort used in all R graphics routines.
#' @param colorC Character. A color, of the sort used in all R graphics routines.
#' @param shapeA Character. Either "c" for circle, "t" for triangle, or "s" for square.
#' @param shapeB Character. See shapeA.
#' @param shapeC Character. See shapeA.
#' @param curves A list of lists of lines. Curves to be plotted.
#' @return NULL.
lineEqualAreaPlotThree <- function(pointsA, pointsB, pointsC, colorA="black", colorB="black", colorC="black", shapeA="c", shapeB="c", shapeC="c", curves=list()) {
  lineEqualAreaPlot(points=c(pointsA, pointsB, pointsC), curves=curves,
                    shapes=c(replicate(length(pointsA), shapeA),
                             replicate(length(pointsB), shapeB),
                             replicate(length(pointsC), shapeC)),
                    colors=c(replicate(length(pointsA), colorA),
                             replicate(length(pointsB), colorB),
                             replicate(length(pointsC), colorC)))
}

#' Equal-area plot of lines, extruded into the third dimension arbitrarily by the user.
#' 
#' @param vs A list of lines.
#' @param zs A vector of real numbers, of the same length as vs. The third coordinate to be plotted.
#' @param ... Other plotting options to be passed to the underlying plot3D.
#' @return NULL.
lineEqualAreaScalarPlot <- function(vs, zs, ...) {
  # Make the points.
  xyzs <- thread(function(v, z) c(equalAreaProjection(lower(v)), z), vs, zs)
  # Make the boundary circle a few times.
  radius <- sqrt(2)
  xys <- lapply((0:72) * 2 * pi / 72, function(theta) {radius * c(cos(theta), sin(theta))})
  first <- lapply(xys, function(xy) c(xy, 0))
  second <- lapply(xys, function(xy) c(xy, min(zs)))
  third <- lapply(xys, function(xy) c(xy, max(zs)))
  plot3D(points=xyzs, curves=list(first, second, third), ...)
}



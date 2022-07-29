


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



#' Equal-area plot of rays.
#'
#' @param points A list of rays, to be plotted as points.
#' @param curves A list of lists of rays. Each list of rays is regarded as a curve to be plotted. Each curve is automatically subidivided at the boundary. Lower-hemisphere parts appear solid, while upper-hemisphere appear dotted.
#' @param colors A character vector. Colors to be passed to the underlying R point-plotting function. They are truncated or recycled as needed.
#' @param shapes A character vector. Shapes to be assigned to the plotted points: "c", "s", or "t" for circle, square, or triangle. Lower-hemisphere (upper-) points are plotted using filled (unfilled) versions of the shapes. They are truncated or recycled as needed.
#' @return NULL.
rayEqualAreaPlot <- function(points=list(), curves=list(), colors=c("black"), shapes=c("c")) {
  # Convert the chosen shapes into their underlying shape codes, considering hemisphere.
  f <- function(i) {
    j <- (i - 1) %% length(shapes) + 1
    if (points[[i]][[3]] > 0) {
      if (shapes[[j]] == "s")
        0
      else if (shapes[[j]] == "t")
        2
      else
        1
    } else {
      if (shapes[[j]] == "s")
        15
      else if (shapes[[j]] == "t")
        17
      else
        19
    }
  }
  underShapes <- sapply(1:length(points), f)
  # All points will be positioned as if on the lower hemisphere.
  underPoints <- lapply(points, lower)
  # Subdivide curves into lower- and upper-hemispherical parts.
  underCurves <- list()
  styles <- c()
  for (curve in curves) {
    curvesSigns <- rayCurvesUpperLower(curve)
    for (i in 1:length(curvesSigns$curves))
      if (curvesSigns$signs[[i]] == 1) {
        underCurves[[length(underCurves) + 1]] <- lapply(curvesSigns$curves[[i]], function(v) -v)
        styles[[length(styles) + 1]] <- "dotted"
      } else {
        underCurves[[length(underCurves) + 1]] <- curvesSigns$curves[[i]]
        styles[[length(styles) + 1]] <- "solid"
      }
  }
  plotEqualArea(points=underPoints, curves=underCurves, colors=colors, shapes=underShapes, styles=styles)
}

#' Equal-area plot of rays, with a circle about each point.
#' 
#' @param points A list of rays, to be plotted as points.
#' @param radii A vector of real numbers. The radii of circles to be drawn about the points. Radii are measured in radians, along the surface of the unit sphere.
#' @param colors A character vector. See rayEqualAreaPlot.
#' @param shapes A character vector. See rayEqualAreaPlot.
#' @return NULL.
rayEqualAreaRadiusPlot <- function(points=list(), radii=c(), colors=c("black"), shapes=c("c")) {
  rayEqualAreaPlot(points, colors=colors, curves=thread(raySmallCircle, points, radii))
}

#' Convenience shortcut for plotting two sets of rays.
#' @param pointsA A list of rays.
#' @param pointsB A list of rays.
#' @param colorA Character. A color, of the sort used in all R graphics routines.
#' @param colorB Character. A color, of the sort used in all R graphics routines.
#' @param shapeA Character. Either "c" for circle, "t" for triangle, or "s" for square.
#' @param shapeB Character. See shapeA.
#' @param curves A list of lists of rays. Curves to be plotted.
#' @return NULL.
rayEqualAreaPlotTwo <- function(pointsA, pointsB, colorA="black", colorB="black", shapeA="c", shapeB="c", curves=list()) {
  rayEqualAreaPlot(points=c(pointsA, pointsB), curves=curves,
                   shapes=c(replicate(length(pointsA), shapeA), replicate(length(pointsB), shapeB)),
                   colors=c(replicate(length(pointsA), colorA), replicate(length(pointsB), colorB)))
}

#' Convenience shortcut for plotting three sets of rays.
#' @param pointsA A list of rays.
#' @param pointsB A list of rays.
#' @param pointsC A list of rays.
#' @param colorA Character. A color, of the sort used in all R graphics routines.
#' @param colorB Character. A color, of the sort used in all R graphics routines.
#' @param colorC Character. A color, of the sort used in all R graphics routines.
#' @param shapeA Character. Either "c" for circle, "t" for triangle, or "s" for square.
#' @param shapeB Character. See shapeA.
#' @param shapeC Character. See shapeA.
#' @param curves A list of lists of rays. Curves to be plotted.
#' @return NULL.
rayEqualAreaPlotThree <- function(pointsA, pointsB, pointsC, colorA="black", colorB="black", colorC="black", shapeA="c", shapeB="c", shapeC="c", curves=list()) {
  rayEqualAreaPlot(points=c(pointsA, pointsB, pointsC), curves=curves,
                   shapes=c(replicate(length(pointsA), shapeA),
                            replicate(length(pointsB), shapeB),
                            replicate(length(pointsC), shapeC)),
                   colors=c(replicate(length(pointsA), colorA),
                            replicate(length(pointsB), colorB),
                            replicate(length(pointsC), colorC)))
}



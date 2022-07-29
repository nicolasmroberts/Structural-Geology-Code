
#Load dependent libraries
library("ggplot2")


# Plot a blank equal area net with some options. 
plotBlankEqualAreaNet <- function(plot = ggplot(), centerCross = TRUE, cardinalTicks = TRUE, cardinalLabels = FALSE, circleLineWeight = 1, fillColor = "white", lineColor = "black") { 
   p = plot
     
   #Draw circle
   border = seq(0,2*pi, length.out = 100)
   circleXY = data.frame("x" = sapply(border, function(s) cos(s)*sqrt(2)), "y" = sapply(border, function(s) sin(s)*sqrt(2)))

   p = plot +
      geom_polygon(data = circleXY, aes(x=x,y=y), color = lineColor, size = circleLineWeight, fill = fillColor) +
      coord_fixed(ratio = 1) +
      theme_void() 

   if(centerCross == TRUE){
      
      center1 = data.frame("x" = c(-0.1, 0.1), "y" = c(0,0))
      center2 = data.frame("x" = c(0,0), "y" = c(-0.1,0.1))
      
      p = p +
         geom_path(data = center1, aes(x=x,y=y), size = circleLineWeight/2, color = lineColor) +   
         geom_path(data = center2, aes(x=x,y=y), size = circleLineWeight/2, color = lineColor)
   }
   
   if(cardinalTicks == TRUE){
      
      tickN = data.frame("x" = c(0, 0), "y" = c(sqrt(2),sqrt(2) + 0.05))
      tickE = data.frame("x" = c(sqrt(2),sqrt(2) + 0.05), "y" = c(0, 0))
      tickS = data.frame("x" = c(0, 0), "y" = c(-sqrt(2), -sqrt(2) - 0.05))
      tickW = data.frame("x" = c(-sqrt(2), -sqrt(2) - 0.05), "y" = c(0, 0))
      
      p = p +
         geom_path(data = tickN, aes(x=x,y=y), size = circleLineWeight, color = lineColor) +
         geom_path(data = tickS, aes(x=x,y=y), size = circleLineWeight, color = lineColor) +
         geom_path(data = tickE, aes(x=x,y=y), size = circleLineWeight, color = lineColor) +
         geom_path(data = tickW, aes(x=x,y=y), size = circleLineWeight, color = lineColor)
      
   }
   
   
   if (cardinalLabels == TRUE){
      
      labels = data.frame(direction = c("N", "E", "S","W"), 
                          x = c(0,sqrt(2) + 0.15, 0, - (sqrt(2) + 0.15)), 
                          y = c(sqrt(2) +0.15, 0, -(sqrt(2) + 0.2), 0)
      )
      
      p = p + 
         geom_text(data = labels, aes(x=x, y=y,label = direction, size = 1, family = "Helvetica"), show.legend = FALSE)
         geom_segment()
   }
   
   return(p)
   
}

plotEqualAreaNetSDTP = function(plot = plotBlankEqualAreaNet(), sdPairs = list(), tpPairs = list(), plotPlanesAsPoles = TRUE, color = "black") {
   p = plot
   
   if(length(sdPairs) > 0){
      if(plotPlanesAsPoles == FALSE) {
         curves = lapply(sdPairs, function(sd) geoCartesianFromStrikeDipDeg(sd))
         p = plotEqualAreaNet(plot = p, curves = list(curves), curveColor = color)
      } else {
         points = lapply(sdPairs, function(sd) geoCartesianFromStrikeDipDeg(sd))
         p = plotEqualAreaNet(plot = p, points = list(points), shape = list(21), color = color)
      }
   }
   
   if(length(tpPairs) > 0) {
      points = lapply(tpPairs, function(tp) geoCartesianFromTrendPlungeDeg(tp))
      p = plotEqualAreaNet(plot = p, points = list(points), color = color)
   }
   return(p)
}


# The main function to call when you want to plot an equal area net with data.
#' @Param plot <- a ggplot object. The default is a blank equal area net, but any equal area net you have already plotted can be called. New data will plot on top of existing data. 
#' @Param points <- a list of lists. Each list is a group of linear features (which plot as points on an EAN)
#' @Param curves <- a list of lists. Each list is a group of planar features (which plot as great circle curves on an EAN)
#' @Param colorBy <- a vector of the same length as your point AND line datasets. Use for gradational color scales, like "northing"
#' @Param edgeColor <- a list of strings with the same length as the points outer list (i.e. number of datasets). Color or the point outlines. 
#' @Param curveColor <- a list of strings. the color of the curves
#' @Param color <- a list of strings. the color fill of the points
#' @Param shape <- a list of numbers correspoinding to the ggplot shapes. The shape of each dataset. 
#' @Param pointSize <- a list of numbers.
#' @Param stroke <- a number. The stroke width of curves
#' @Param pointAlpha <- a list of numbers between 0 and 1
#' @Param curveAlpha <- a list of numbers between 0 and 1
#' @Param numCurveSegs <- a number. The number of straight segments that are used to draw each great circle curve
#' @Param showLegend
plotEqualAreaNet = function(plot = plotBlankEqualAreaNet(), points = NULL, curves = NULL, colorBy = NULL, edgeColor = list("white"), curveColor = list("black"), color = list("black"), shape = list(22), pointSize = list(3), stroke = 0.5, pointAlpha = list(1), curveAlpha = list(1), numCurveSegs = 80,  showLegend = FALSE){

   p = plot

   
   if (is.null(curves) == FALSE) {
      if (length(curveColor) != length(curves)) {
         curveColor = rep(curveColor[[1]], length(curves))
      }
      if (length(curveAlpha) != length(curves)) {
         curveAlpha = rep(curveAlpha[[1]], length(curves))
      }
      for (i in 1:length(curves)) {
         curvesXY = curvesXYFromPole(curves[[i]], numCurveSegs = numCurveSegs)
         if (is.null(colorBy) == FALSE){
            colorByCurve = rep(colorBy, each = nrow(curvesXY[curvesXY$index == 1,]))
            p = p +
               geom_path(data = curvesXY,  aes(x = X1, y = X2, group = index, color = colorByCurve), size = stroke, alpha = curveAlpha[[i]], show.legend = showLegend) 
         } else {
            p = p +
               geom_path(data = curvesXY,  aes(x = X1, y = X2, group = index), color = curveColor[[i]], size = stroke, alpha = curveAlpha[[i]], show.legend = showLegend) 
         }
      }
      
   }
   
   
   
   if (is.null(points) == FALSE) {
      if (length(shape) != length(points)) {
         shape = rep(shape[[1]], length(points))
      }
      if (length(color) != length(points)) {
         color = rep(color[[1]], length(points))
      }
      if (length(pointAlpha) != length(points)) {
         pointAlpha = rep(pointAlpha[[1]], length(points))
      }
      if (length(pointSize) != length(points)) {
         pointSize = rep(pointSize[[1]], length(points))
      }
      if (length(edgeColor) != length(points)) {
         edgeColor = rep(edgeColor[[1]], length(points))
      }
      for (i in 1:length(points)){
         linesXY = data.frame(t(sapply(points[[i]], function(s) equalAreaProjection(lower(s)))))
         names(linesXY) = c("linX", "linY")
            if (is.null(colorBy) == FALSE){
               p = p + 
               geom_point(data = linesXY, aes(x= linX, y=linY, fill = colorBy), color = edgeColor[[i]], shape = shape[[i]], size = pointSize[[i]], alpha = pointAlpha[[i]], show.legend = showLegend)
               }else{
               p = p + 
               geom_point(data = linesXY, aes(x= linX, y=linY), fill = color[[i]], color = edgeColor[[i]], shape = shape[[i]], size = pointSize[[i]], alpha = pointAlpha[[i]], show.legend = showLegend)
            }
         }
   }
p
}

plotEqualAreaCI = function(plot = plotBlankEqualAreaNet(), polygonPoints = list(), edgeColor = list("black"), edgeSize = 1, fillColor = list("grey75"), alpha = 1) {
   p = plot
   
   if(length(edgeColor) < length(polygonPoints)){
      edgeColor = rep(edgeColor[[1]], length(polygonPoints))
   }
   if(length(fillColor) < length(polygonPoints)){
      fillColor = rep(fillColor[[1]], length(polygonPoints))
   }
    
   for(i in 1:length(polygonPoints)){
      linesXY = data.frame(t(sapply(polygonPoints[[i]], function(s) equalAreaProjection(lower(s)))))
      names(linesXY) = c("linX", "linY")
      p = p +
         geom_polygon(data = linesXY, aes(x = linX, y = linY), fill = fillColor[[i]], size = edgeSize, color = edgeColor[[i]], alpha = alpha)
   }  
   p
}

curvesXYFromPole = function(poles, numCurveSegs = 100) {
   curves = lapply(poles, function(p) rayGreatCircleHorizontalStart(p, numCurveSegs))
   
   for(i in 1:length(curves)){
      curveTemp = data.frame(t(sapply(curves[[i]], function(s) s)))
      #curves[[i]] = curveTemp[curveTemp$X3 < 0,]
      curves[[i]] = curveTemp
   }
   
   
   curvesXY = lapply(curves, function(p) data.frame(t(apply(p,1, function(g) equalAreaProjection(lower(g))))))
   
   curvesXYDataFrame = data.frame(t(sapply(1:nrow(curvesXY[[1]]), function(s) curvesXY[[1]][s,])))
   curvesXYDataFrame$index = 1
   
   for(i in 2:(length(curvesXY))) {
      
      tempFrame = data.frame(t(sapply(1:nrow(curvesXY[[i]]), function(s) curvesXY[[i]][s,])))
      tempFrame$index = i
      
      curvesXYDataFrame = rbind(curvesXYDataFrame, tempFrame)
   }
   
   curvesXYDataFrame$X1 = unlist(curvesXYDataFrame$X1)
   curvesXYDataFrame$X2 = unlist(curvesXYDataFrame$X2)
   curvesXYDataFrame
}




#' A ray perpendicular to a given vector, deterministically.
#' 
#' The result is deterministic but arbitrary. Due to the hairy ball theorem, the result cannot be chosen to depend smoothly on the input.
#' @param v A 3-dimensional vector. Need not be unit.
#' @return A ray, perpendicular to v.
rayOrthogonalHorizontal <- function(v) {
     strikeDip = geoStrikeDipDegFromCartesian(v)
     horizontal <- geoCartesianFromTrendPlungeDeg(c(strikeDip[1], 0))
     horizontal
}

#' The great circle perpendicular to the given ray.
#' 
#' @param pole A ray.
#' @param numSteps The number of line segments to be used in the approximation of the great circle.
#' @return A list of numSteps+1 rays. The first and last one are identical. Otherwise, they are evenly spaced on the circle perpendicular to the given pole, and in order (either clockwise or counter-clockwise, depending on your viewpoint).
rayGreatCircleHorizontalStart <- function(pole, numSteps=50) {
     v <- rayOrthogonalHorizontal(pole)
     w <- cross(pole, v)
     angles <- (1:numSteps) * (-pi / (numSteps + 1) )
     angles <- c(-10^(-15), angles, -pi + 10^-15)
     lapply(angles, function(a) {cos(a) * v + sin(a) * w})
}





lineFunctionPlotBoxSets <- function(f, levels, numNonAdapt, ...) {
     # Non-adaptive refinement.
     boxes <- lineFunctionPlotBase(f, ...)
     if (numNonAdapt >= 1)
          for (i in 1:numNonAdapt) 
               boxes <- unlist(lapply(boxes, lineFunctionPlotRefinement, f, ...), recursive=FALSE)
     boxes
}

# lineFunctionPlotBoxSets customized to the problem of Kamb contouring.
lineKambBoxSets <- function(points, multiples=c(3, 6, 9, 12), k=3, degree=3, numNonAdapt=4) {
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
     lineFunctionPlotBoxSets(f, levels, numNonAdapt)
}




#Plot an equal area net with kamb contours (no point data will be plotted. To plot point data on top of kamb contours, save this plot as a variable and use it as the plot in the plotEqualAreaNet function)
plotKambContours = function(kambPlot = ggplot(), points = list(), breaks = c(0,3,6,9,12), numNonAdapt = 4, showLegend = FALSE){
     
     boxes = lineKambBoxSets(points, numNonAdapt = numNonAdapt, multiples = breaks)
     density = c()
     
     for(i in 1: length(boxes)){
          density[i] = mean(sapply(boxes[[i]], function(s) s[4]))
     }
     
     breaks = c(breaks, max(density))
     
     boxXYs = data.frame(t(sapply(boxes[[1]], function(s) equalAreaProjection(s[1:3]))))
     boxXYs$density = cut(mean(sapply(boxes[[1]], function(s) s[4])), breaks = breaks, include.lowest = TRUE)
     boxXYs$id = 1
     
     for(i in 2:length(boxes)) {
          boxXY = data.frame(t(sapply(boxes[[i]], function(s) equalAreaProjection(s[1:3]))))
          boxXY$density = cut(mean(sapply(boxes[[i]], function(s) s[4])), breaks = breaks, include.lowest = TRUE)
          boxXY$id = i
          boxXYs = rbind(boxXYs, boxXY)
     }
     
     kambPlot = kambPlot + geom_polygon(data = boxXYs, aes(x=X1, y = X2, fill = factor(density), color = factor(density), group = id), show.legend = showLegend)
     
     kambPlot = kambPlot + scale_fill_brewer(palette = 6) + scale_color_brewer(palette = 6) 
     
     
     if(showLegend == TRUE){
          kambPlot = kambPlot +
               theme(legend.position = "right") + labs(fill = expression(sigma), color = expression(sigma)) + guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) 
     }
     kambPlot = plotBlankEqualAreaNet(plot = kambPlot, fillColor = alpha("white", 0))
     kambPlot
}






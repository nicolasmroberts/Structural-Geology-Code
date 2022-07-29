rotateSymbol = function(sdPolygon, angleDeg, halfSize = 1){
     theta = angleDeg*degree
     #sdLine = sdLine*halfSize
     sdPolygon = sdPolygon*halfSize
     rotation = rbind(c(cos(theta), sin(theta)),
                      c(-sin(theta), cos(theta)))
     # sdLineNew = data.frame(t(sapply(1:nrow(sdLine), function(i) rotation%*%t(sdLine[i,]))))
     sdPolygonNew = data.frame(t(sapply(1:nrow(sdPolygon), function(i) rotation%*%t(sdPolygon[i,]))))
     return(sdPolygonNew)
}


plotMapStructureSymbols = function(plot = ggplot(), lat = NULL, long = NULL, strikeDips = NULL, trendPlunges = NULL, lineCartesian = NULL, poleCartesian = NULL, colorBy = NULL, pointSize = 4, symbolColor = "black") {
     # Set the size of the symbols based on the map span. 
     xSpan = max(long) - min(long)
     
     
     # Define symbols
     sdTriangle = data.frame(x = c(0,0,0,.25,0,0),
                             y = c(-1,1,0.25,0,-.25,-1))
     
     tpPolygon = data.frame(x =c(0,0,.2,-.2,0,0), 
                            y =c(-1,1.4,1,1,1.4,0))
     
     # Get xy coordinates into a dataframe
     dataFrame = data.frame(cbind(long, lat))
     
     # Append strike-dip information to dataFrame
     if(is.null(strikeDips) == FALSE){
          dataFrame$sd = strikeDips
     }
     if(is.null(poleCartesian) == FALSE){
          dataFrame$sd = lapply(poleCartesian, function(vec) geoStrikeDipDegFromCartesian(vec))
     }
     if(is.null(trendPlunges) == FALSE){
          dataFrame$tp = trendPlunges
     }
     if(is.null(lineCartesian) == FALSE){
          dataFrame$tp = lapply(lineCartesian, function(vec) geoTrendPlungeDegFromCartesian(vec))
     }
     
     
     #If coloring by a parameter, plot colored circles on the plot.
     if (is.null(colorBy) == FALSE) {
          if(colorBy == "dip"){
               Dips = sapply(dataFrame$sd, function(sd) sd[2])
               plot = plot + geom_point(data = dataFrame, aes(x = long, y = lat, color = Dips), size = pointSize)
          }
          else if(colorBy == "plunge"){
               Plunges = sapply(dataFrame$tp, function(tp) tp[2])
               plot = plot + geom_point(data = dataFrame, aes(x = long, y = lat, color = Plunges), size = pointSize)
          } else {
               plot = plot + geom_point(data = dataFrame, aes(x = long, y = lat, color = colorBy), size = pointSize)
          }
          plot = plot + scale_color_viridis(limits = c(0,90), direction = -1)
     }
     
     #Finally, plot the strike-dip and/or trend-plunge symbols. 
     if (is.null(poleCartesian) == FALSE | is.null(strikeDips) == FALSE) {
          #rotate the symbols and plot them, one by one
          for(i in 1:length(dataFrame$sd)){
               sd = dataFrame$sd[[i]]
               symbol = rotateSymbol(sdTriangle, sd[1], halfSize = xSpan/30)
               # symbol[[1]]$X1 = symbol[[1]]$X1 + xyCoords$x[i]
               # symbol[[1]]$X2 = symbol[[1]]$X2 + xyCoords$y[i]
               symbol$X1 = symbol$X1 + dataFrame$long[i]
               symbol$X2 = symbol$X2 + dataFrame$lat[i]
               plot = plot + 
                    #geom_line(data = symbol[[1]],aes(x = X1, y =X2), color = symbolColor) +
                    geom_polygon(data = symbol, aes(x=X1, y=X2), color = symbolColor, fill = symbolColor, size = .5)
          }
          plot = plot + 
               #theme_void() + 
               theme(panel.grid.major = element_blank(), 
                     legend.position = c(1,0),
                     legend.margin = margin(0.005,.005,.005,.005, "npc"),
                     legend.justification = c(1,0),
                     legend.key.height = unit(0.04, "npc"), 
                     legend.direction = "horizontal",
                     legend.background = element_rect(fill = "white", color = "black", size = .5),
                     panel.border = element_rect(size = 1, color = "black", fill = NA))
     }
     if (is.null(lineCartesian) == FALSE | is.null(trendPlunges) == FALSE) {
          
          
          #rotate the symbols and plot them, one by one
          for(i in 1:length(dataFrame$tp)){
               tp = dataFrame$tp[[i]]
               arrowAdj = (tp[2]/90)
               tpPolygon = data.frame(x =c(0,0,.2,-.2,0,0), 
                                      y =c(-1 + arrowAdj,1.4-arrowAdj,1-arrowAdj,1 -arrowAdj,1.4 - arrowAdj,-1 + arrowAdj))
               symbol = rotateSymbol(tpPolygon, tp[1], halfSize = xSpan/30)
               symbol$X1 = symbol$X1 + dataFrame$long[i]
               symbol$X2 = symbol$X2 + dataFrame$lat[i]
               plot = plot + 
                    #geom_line(data = symbol[[1]],aes(x = X1, y =X2), color = symbolColor) +
                    geom_polygon(data = symbol, aes(x=X1, y=X2), color = symbolColor, fill = symbolColor)
          }
          
          plot = plot + 
               #theme_void() + 
               theme(panel.grid.major = element_blank(), 
                     legend.position = c(1,0),
                     legend.margin = margin(0.005,.005,.005,.005, "npc"),
                     legend.justification = c(1,0),
                     legend.key.height = unit(0.04, "npc"), 
                     legend.direction = "horizontal",
                     legend.background = element_rect(fill = "white", color = "black", size = .5),
                     panel.border = element_rect(size = 1, color = "black", fill = NA))
     }  
     
     
     plot + coord_fixed()
     
     
     
     
}


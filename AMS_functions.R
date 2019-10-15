


#this loads a FUNCTION that appends an ams datset with easting and northing information from the locData variable.
#inputs: it takes an ams dataset, a list of station names with easting and northings, and the number of letters/numbers in each station name. For example, if the station name is AME16_037, then the number would be 9. Note that all your station names should be the same length. so station nubmers should be 003, 025, and 100 rather than 3, 25, and 100.
#' locationPairing
#'
#' @param amsData. A data frame that was output by the geoEllipsoidDataFromAGICOFile function in Josh Davis' geologyGeometry library.
#' @param locData. A dataframe that has columns labelled "name", "easting", and "northing".
#' @param uniqueStation. Scalar quantity. The number of characters in each station name. 
#'
#' @return amsData with easting and northing columns tacked on. 

locationPairing <- function(amsData, locData, uniqueStation){
     amsData$easting <- replicate(nrow(amsData), 1)
     amsData$northing <- replicate(nrow(amsData), 1)
     amsData$station <- replicate(nrow(amsData), "")
     for(i in 1:nrow(amsData)){
          name <- substr(amsData$Name[i], 1, uniqueStation)
          locDataSub = subset(locData, grepl(paste(name), locData$station))
          amsData$easting[i] <- locDataSub$easting[1]
          amsData$northing[i] <- locDataSub$northing[1]
          amsData$station[i] <- name
     }
     return(amsData)  
}


stationDataSubset = function(stationName, amsData) {
     stationData <- subset(amsData, grepl(paste(stationName), amsData$station))
     return(stationData)
}

getKsFromAmsData <- function(stationData, whichAxes = "all") {
     if (whichAxes == "all") {
     K1 = data.frame(t(sapply(stationData$rotation, function(s) equalAreaProjection(lower(s[1,])))))
     K2 = data.frame(t(sapply(stationData$rotation, function(s) equalAreaProjection(lower(s[2,])))))
     K3 = data.frame(t(sapply(stationData$rotation, function(s) equalAreaProjection(lower(s[3,])))))
     
     xy = data.frame(rbind(K1,K2,K3))
     xy = data.frame(xy, "axis"=rep(c(1,2,3), each=length(stationData$Name)))
     
     }
     
     if (whichAxes == "K1") {
          K1 = data.frame(t(sapply(stationData$rotation, function(s) equalAreaProjection(lower(s[1,])))))
          xy = data.frame(rbind(K1))
          xy = data.frame(xy, "axis"=rep(c(1), length(stationData$Name)))
     }
     
     if (whichAxes == "K3") {
          K3 = data.frame(t(sapply(stationData$rotation, function(s) equalAreaProjection(lower(s[3,])))))
          xy = data.frame(rbind(K3))
          xy = data.frame(xy, "axis"=rep(c(3), length(stationData$Name)))
     }
     return(xy)
}

getKsFromAmsMean <- function(dataMean) {
     Kmean1 = data.frame(t(equalAreaProjection(lower(dataMean$rotation[3,]))))
     Kmean2 = data.frame(t(equalAreaProjection(lower(dataMean$rotation[2,]))))
     Kmean3 = data.frame(t(equalAreaProjection(lower(dataMean$rotation[1,]))))
     
     meanXY = data.frame(rbind(Kmean1, Kmean2, Kmean3))
     meanXY = data.frame(meanXY, "axis" = rep(c(1,2,3), each = 1))
     
     return(meanXY)
}


plotNet <- function(amsData, showFollin = TRUE, whichAxes = "all") { 
     
     stationData <- amsData
     xy <- getKsFromAmsData(stationData, whichAxes)
     
     dataMean <- ellMean(stationData$vector)
     
     plane=rayGreatCircle(dataMean$rotation[1,], 500)
     planeXY = data.frame(t(sapply(plane, function(s) equalAreaProjection(lower(s)))))
     lineXY = data.frame(t(equalAreaProjection(lower(dataMean$rotation[3,]))))
     
     
     shapeVector = c("1" = 22 , "2" = 24 , "3" = 21)
     colorVector = c("1" = "blue", "2" = "green", "3" = "#e75480")
     n = length(stationData$Name)
     
     return(
          
     if (showFollin == TRUE) {
          ggplot() +
               geom_point(data=planeXY, aes(x=X1, y=X2), color = "black", size = 1, show.legend = FALSE) +
               geom_point(data=lineXY, aes(x=X1, y=X2), color = "black", size = 5, show.legend = FALSE) +
               geom_point(data=xy, aes(x=X1, y=X2, shape = as.factor(axis), fill=as.factor(axis)), size = 2, show.legend = FALSE) +
               geom_circle(aes(x0=0,y0=0,r=sqrt(2)), size = 1, inherit.aes = FALSE) +
               coord_fixed(ratio = 1) +
               scale_shape_manual(values = shapeVector) +
               scale_fill_manual(values = colorVector) +
               theme_void() +
               labs(title = paste0(stationData$station, " (n = ", n,")"))
     } else {
          ggplot() +
               geom_point(data=xy, aes(x=X1, y=X2, shape = as.factor(axis), fill=as.factor(axis)), size = 2, show.legend = FALSE) +
               geom_circle(aes(x0=0,y0=0,r=sqrt(2)), size = 1, inherit.aes = FALSE) +
               coord_fixed(ratio = 1) +
               scale_shape_manual(values = shapeVector) +
               scale_fill_manual(values = colorVector) +
               theme_void() +
               labs(title = paste0(stationData$station, " (n = ", n,")"))
     }
     )
}

plotNetStation <- function(stationName, amsData) { 
     
     stationData <- stationDataSubset(stationName, amsData)
     plotNet(stationData)
}


plotHsuStation = function(stationName, amsData) {
     stationData <- stationDataSubset(stationName, amsData)
     return( plotHsu(stationData))
}


plotHsu = function(amsData) {
stationData <- amsData
stationMean <- amsMean(amsData)
logAs <- lapply(stationData$logA, function(ell) ell)



x <- function(logA) {-sqrt(1.5) * max(logA - sum(logA) / 3) - sqrt(1.5) * min(logA - sum(logA) / 3)}
y <- function(logA) {sqrt(0.5) * max(logA - sum(logA) / 3) - sqrt(0.5) * min(logA - sum(logA) / 3)}

lodes <- sapply(logAs, x)
oct <- sapply(logAs, y)
meanLode <- sapply(stationMean$logA, x)
meanOct <- sapply(stationMean$logA, y)

# the maximum octohedral shear strain
maxEs <- max(oct + 0.05)

# Dataframe that stores the Hsu Plot line segments that show lines of equal Lodes parameter
lines <- data.frame(x1 = c(0,0,0), 
                    x2 = c(maxEs*cos(pi/3),maxEs*cos(pi/2), maxEs*cos(2/3 * pi)),
                    y1 = c(0,0,0), 
                    y2 = c(maxEs*sin(pi/3), maxEs*sin(pi/2), maxEs*sin(2/3 * pi)),
                    label = c(1, 0, -1),
                    lineAlpha = c(1,.75,1)
)

# Make data frame of Hsu plot arcs of equal oct. shear strain
# Change arcInterval to have more or fewer arcs display on the Hsu plot
arcInterval = 0.1

arcs <- data.frame(
     start = rep(-pi/6, maxEs/arcInterval + 1),
     end = rep(pi/6, maxEs/arcInterval + 1),
     r = seq(from = 0, to = maxEs , by = arcInterval),
     textX = seq(from = 0*cos(pi/3), to = maxEs*cos(pi/3), by=arcInterval * cos(pi/3)),
     textY = seq(from = 0*sin(pi/3), to = maxEs*sin(pi/3), by=arcInterval * sin(pi/3))
)

#construct dataframe to plot
xy <- data.frame(cbind(oct, lodes))
xyMean <- data.frame(cbind(meanOct, meanLode))

#The Plot
HsuPlot <-ggplot() +
     coord_fixed(ratio = 1, ylim = c(-0.01,maxEs + 0.02), xlim = c(-maxEs/2 -0.01,maxEs/2 + 0.01)) +
     #plot the arcs and label them
     geom_arc(aes(x0 = 0, y0 = 0, r = r, start = start, end = end, alpha = 0.9), linetype = 1, color = "dark grey", data = arcs, show.legend = FALSE) +
     geom_text(aes(x=textX, y=textY,label = r,angle = -30), size = 2.75, family = "Helvetica", nudge_x = 0.01, nudge_y = -0.01, data = arcs)+
     geom_arc(aes(x0=0,y0=0, r= maxEs, start = - pi/6, end = pi/6), linetype = 1, show.legend = FALSE) +     
     #plot lines of equal lodes parameter and label them
     geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, alpha = lineAlpha), linetype = 1, data = lines, show.legend = FALSE) +
     geom_text(aes(x=x2, y=y2,label = label, vjust = -0.5, family = "Helvetica"), size = 2.75, data = lines) +
     #geom_text(aes(x=x, y=y, label = labels, angle = angle, family = "Times"), size = 2.75, data = axisLabels)+
     #plot data
     geom_point(aes(x=lodes, y=oct), data = xy, size = 2)+
     geom_point(aes(x=meanLode, y=meanOct), data = xyMean, size = 4, color = "blue")+
     theme_void() +
     theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

return(HsuPlot)
}



amsMean = function(amsData) {
     amsStationData <-amsData
     meanEll <- ellMean(amsStationData$vector)
     rot <- meanEll$rotation
     axis <- meanEll$a
     logAxis <- meanEll$logA
     #Store the mean data for each iteration
     
     stationMean <- data.frame(
          Name = ("MeanData"),
          easting = (as.numeric(amsStationData$easting[[1]])),
          northing = (as.numeric(amsStationData$northing[[1]])),
          specimens = (nrow(amsStationData)),
          vector = t(list(meanEll$vector)),
          # # sorting of the axis lengths/rotation so that it matches amsData
          logA = t(list(sort(logAxis, decreasing=TRUE))),
          a = t(list(sort(axis, decreasing=TRUE))),
          rotation = t(list(rbind(rot[3,], rot[2,], rot[1,]))),
          tensor = t(list(meanEll$tensor)),
          # mean susceptibility for the station
          Km = (mean(unlist(amsStationData$Km)))
     )
     
     return(stationMean)
     # #return(amsData)

}


stationMean = function(stationName, amsData) {
     amsStationData <- subset(amsData, grepl(paste(stationName), amsData$station))
     return(amsMean(amsStationData))
}



plotNetStationMean <- function(stationMean) { 
     
     xy <- getKsFromAmsMean(stationMean)
     
     shapeVector = c("1" = 22 , "2" = 24 , "3" = 21)
     colorVector = c("1" = "blue", "2" = "green", "3" = "#e75480")
     
     n = length(stationData$Name)
     
     
     return(
          
          ggplot() +
               geom_point(data=xy, aes(x=X1, y=X2, shape = as.factor(axis), fill=as.factor(axis)), color = "black", size = 4, show.legend = FALSE) +
               geom_circle(aes(x0=0,y0=0,r=sqrt(2)), inherit.aes = FALSE) +
               coord_fixed(ratio = 1) +
               scale_shape_manual(values = shapeVector) +
               scale_color_manual(values = colorVector) +
               theme_void() +
               labs(title = paste0(stationData$station, " (n = ", n,")"))
          
     )
     
}


computeAMSBootStrap95 <- function(amsData) {
     data <- amsData
     
     if (nrow(data) > 5) {
          boots <- ellBootstrapInference(data$vector, 1000)
          boots95 <- ellHighEllipsoidVectors(boots$covarInv, boots$center, (boots$q095)^2)
          
          bootsEllFrame <- data.frame()
          
          for (i in 1:length(boots$bootstraps)) {
               bootEll <- ellEllipsoidFromVector(boots$bootstraps[[i]])
               bootEllFrame <- data.frame()
               frame <- data.frame(
                    Name = ("Bootstraps"),
                    vector = t(list(bootEll$vector)),
                    # # sorting of the axis lengths/rotation so that it matches amsData
                    logA = t(list(sort(bootEll$logA, decreasing=TRUE))),
                    a = t(list(sort(bootEll$a, decreasing=TRUE))),
                    rotation = t(list(rbind(bootEll$rotation[3,], bootEll$rotation[2,], bootEll$rotation[1,]))),
                    tensor = t(list(bootEll$tensor))
               )
               
               bootsEllFrame <- rbind(bootsEllFrame, frame) 
          }
          return(bootsEllFrame)
     }
     else{
          return("not enough data points")
     }
}




computeStationAMSBootStrap95 <- function(stationName, amsData) {
     data <- stationDataSubset(stationName, amsData)
     computeAMSBootStrap95(data)
     
}


plotNetBootsAll <- function(amsData, boots = NULL, showFollin = TRUE) {
     if (is.null(boots) == TRUE) {
          boots <- computeAMSBootStrap95(amsData)
     }
data <- amsData
dataMean <- ellMean(data$vector)

plane=rayGreatCircle(dataMean$rotation[1,], 500)
planeXY = data.frame(t(sapply(plane, function(s) equalAreaProjection(lower(s)))))
# lineXY = data.frame(t(equalAreaProjection(lower(dataMean$rotation[3,]))))

xy = getKsFromAmsData(data)
xyBoots = getKsFromAmsData(boots)
meanXY = getKsFromAmsMean(dataMean)

shapeVector = c("1" = 22 , "2" = 24 , "3" = 21)
colorVector = c("1" = "blue", "2" = "green", "3" = "#e75480")

n = length(data$Name)

return(
     
     if (showFollin == TRUE) {
     ggplot() +
          geom_point(data=planeXY, aes(x=X1, y=X2), color = "black", size = 1, show.legend = FALSE) +
          geom_point(data=xyBoots, aes(x=X1, y=X2, shape = as.factor(axis), fill=as.factor(axis)), stroke = 0, size = 1, alpha = 0.5, show.legend = FALSE) +
          geom_point(data=xy, aes(x=X1, y=X2, shape = as.factor(axis), fill=as.factor(axis)), stroke = 1, size = 3, show.legend = FALSE) +
          # geom_point(data=lineXY, aes(x=X1, y=X2), shape = 22, color = "black", fill = "white", size = 5, show.legend = FALSE) +
          # geom_point(data=meanXY, aes(x=X1, y=X2, shape = as.factor(axis)), fill="black", color = "black", size = 5, show.legend = FALSE) +
          geom_circle(aes(x0=0,y0=0,r=sqrt(2)), size = 1, inherit.aes = FALSE) +
          coord_fixed(ratio = 1) +
          scale_shape_manual(values = shapeVector) +
          scale_fill_manual(values = colorVector) +
          theme_void() +
          labs(title = paste0(boots$Name[[1]], ", n = ",n), element_text(size = 10, family = "Helvetica"))
     }
     else {
          geom_point(data=planeXY, aes(x=X1, y=X2), color = "black", size = 1, show.legend = FALSE) +
               geom_point(data=xyBoots, aes(x=X1, y=X2, shape = as.factor(axis), fill=as.factor(axis)), stroke = 0, size = 1, alpha = 0.5, show.legend = FALSE) +
               geom_point(data=xy, aes(x=X1, y=X2, shape = as.factor(axis), fill=as.factor(axis)), stroke = 1, size = 3, show.legend = FALSE) +
               # geom_point(data=lineXY, aes(x=X1, y=X2), shape = 22, color = "black", fill = "white", size = 5, show.legend = FALSE) +
               #geom_point(data=meanXY, aes(x=X1, y=X2, shape = as.factor(axis)), fill="black", color = "black", size = 5, show.legend = FALSE) +
               geom_circle(aes(x0=0,y0=0,r=sqrt(2)), size = 1, inherit.aes = FALSE) +
               coord_fixed(ratio = 1) +
               scale_shape_manual(values = shapeVector) +
               scale_fill_manual(values = colorVector) +
               theme_void() +
               labs(title = paste0(boots$Name[[1]], ", n = ",n), element_text(size = 10, family = "Helvetica"))
     } )
}


plotNetBoots <- function(amsData, stationName, boots = NULL) { 
     if (is.null(boots) == TRUE) {
          boots <- computeAMSBootStrap95(stationName, amsData)
     }
     data <- stationDataSubset(stationName, amsData)
     dataMean <- ellMean(data$vector)
     
     plane=rayGreatCircle(dataMean$rotation[1,], 500)
     planeXY = data.frame(t(sapply(plane, function(s) equalAreaProjection(lower(s)))))
     # lineXY = data.frame(t(equalAreaProjection(lower(dataMean$rotation[3,]))))
     
     
     xy = getKsFromAmsData(data)
     xyBoots = getKsFromAmsData(boots)
     meanXY = getKsFromAmsMean(dataMean)
     
     shapeVector = c("1" = 22 , "2" = 24 , "3" = 21)
     colorVector = c("1" = "blue", "2" = "green", "3" = "#e75480")
     
     n = length(data$Name)
     
     return(
          
          ggplot() +
               geom_point(data=planeXY, aes(x=X1, y=X2), color = "black", size = 1, show.legend = FALSE) +
               geom_point(data=xyBoots, aes(x=X1, y=X2, shape = as.factor(axis), fill=as.factor(axis)), stroke = 0, size = 1, alpha = 0.5, show.legend = FALSE) +
               geom_point(data=xy, aes(x=X1, y=X2, shape = as.factor(axis), fill=as.factor(axis)), stroke = 1, size = 3, show.legend = FALSE) +
               # geom_point(data=lineXY, aes(x=X1, y=X2), shape = 22, color = "black", fill = "white", size = 5, show.legend = FALSE) +
               #geom_point(data=meanXY, aes(x=X1, y=X2, shape = as.factor(axis)), fill="black", color = "black", size = 5, show.legend = FALSE) +
               geom_circle(aes(x0=0,y0=0,r=sqrt(2)), size = 1, inherit.aes = FALSE) +
               coord_fixed(ratio = 1) +
               scale_shape_manual(values = shapeVector) +
               scale_fill_manual(values = colorVector) +
               theme_void() +
               labs(title = paste0(boots$Name[[1]], ", n = ",n), element_text(size = 10, family = "Helvetica"))
      
          
          
     )
}


computeAndPlotAllBootstraps <- function(amsData, stationList, path, fileType = "pdf") {
     locationNames <- list()
     bootPlots <- list()
     hsuPlots <- list()
     
     for (i in 1:length(stationList)){
          print(paste(i, " of ", length(stationList[[i]])))
          locationNames[[i]] = stationList[[i]]
          subset = stationDataSubset(stationList[[i]], amsData)
          if (length(subset$Name) <= 5) {
               
               bootPlots[[i]] <- plotNetStation(stationList[[i]], amsData)
               hsuPlots[[i]] <- plotHsuStation(stationList[[i]], amsData)
               
          }
          else {
               boots <- computeAMSBootStrap95(stationList[[i]], amsData)
               bootPlots[[i]] <- plotNetBoots(amsData, stationList[[i]], boots)
               hsuPlots[[i]] <- plotBootstrapsHsu(stationList[[i]], amsData, boots)
          }
     }
     plotsAndNames <- cbind(bootPlots, hsuPlots, locationNames)
     for (i in 1:length(locationNames)){
          plots <- list(plotsAndNames[i,]$bootPlots, plotsAndNames[i,]$hsuPlots)
          lay = rbind(c(1,2))
          
          if (fileType == "png") {
          ggsave(paste0(path, "/", plotsAndNames[i,]$locationNames, ".png"), plot = marrangeGrob(grobs = plots, layout_matrix = lay), height = 3, width = 6, unit = "in", dpi = 60)
          }
          else {
               ggsave(paste0(path, "/", plotsAndNames[i,]$locationNames, ".", fileType), plot = marrangeGrob(grobs = plots, layout_matrix = lay), height = 4, width = 8, unit = "in")
               
          }
          
          }
     
     
     
     return(plotsAndNames)
}


plotBootstrapsHsu = function(stationName,amsData, boots = NULL) {
     if (is.null(boots) == TRUE) {
          boots <- computeAMSBootStrap95(stationName, amsData)
     }
     stationData <- stationDataSubset(stationName, amsData)
     stationMean <- stationMean(stationName, amsData)
     
     logAs <- lapply(stationData$logA, function(ell) ell)
     
     bootStraps <- boots$vector
     
     getLogA = function(vec) {
        ellipsoid = ellEllipsoidFromVector(vec)
        return(ellipsoid$logA)
     } 
     logAsBoots <- lapply(bootStraps, getLogA)
     
     x <- function(logA) {-sqrt(1.5) * max(logA - sum(logA) / 3) - sqrt(1.5) * min(logA - sum(logA) / 3)}
     y <- function(logA) {sqrt(0.5) * max(logA - sum(logA) / 3) - sqrt(0.5) * min(logA - sum(logA) / 3)}
     
     lodes <- sapply(logAs, x)
     oct <- sapply(logAs, y)
     lodesBoots <- sapply(logAsBoots, x)
     octBoots <- sapply(logAsBoots, y)
     meanLode <- sapply(stationMean$logA, x)
     meanOct <- sapply(stationMean$logA, y)
     
     # the maximum octohedral shear strain
     maxEs <- max(oct + 0.05)
     
     # Dataframe that stores the Hsu Plot line segments that show lines of equal Lodes parameter
     lines <- data.frame(x1 = c(0,0,0), 
                         x2 = c(maxEs*cos(pi/3),maxEs*cos(pi/2), maxEs*cos(2/3 * pi)),
                         y1 = c(0,0,0), 
                         y2 = c(maxEs*sin(pi/3), maxEs*sin(pi/2), maxEs*sin(2/3 * pi)),
                         label = c(1, 0, -1),
                         lineAlpha = c(1,.75,1)
     )
     
     # Make data frame of Hsu plot arcs of equal oct. shear strain
     # Change arcInterval to have more or fewer arcs display on the Hsu plot
     arcInterval = 0.1
     
     arcs <- data.frame(
          start = rep(-pi/6, maxEs/arcInterval + 1),
          end = rep(pi/6, maxEs/arcInterval + 1),
          r = seq(from = 0, to = maxEs , by = arcInterval),
          textX = seq(from = 0*cos(pi/3), to = maxEs*cos(pi/3), by=arcInterval * cos(pi/3)),
          textY = seq(from = 0*sin(pi/3), to = maxEs*sin(pi/3), by=arcInterval * sin(pi/3))
     )
     
     #construct dataframe to plot
     xy <- data.frame(cbind(oct, lodes))
     xyMean <- data.frame(cbind(meanOct, meanLode))
     xyBoots <- data.frame(cbind(octBoots, lodesBoots))
     
     #The Plot
     HsuPlot <- ggplot() +
          coord_fixed(ratio = 1, ylim = c(-0.01,maxEs + 0.02), xlim = c(-maxEs/2 -0.01,maxEs/2 + 0.01)) +
          #plot the arcs and label them
          geom_arc(aes(x0 = 0, y0 = 0, r = r, start = start, end = end, alpha = 0.9), linetype = 1, color = "dark grey", data = arcs, show.legend = FALSE) +
          geom_text(aes(x=textX, y=textY,label = r,angle = -30), size = 2.75, family = "Helvetica", nudge_x = 0.01, nudge_y = -0.01, data = arcs)+
          geom_arc(aes(x0=0,y0=0, r= maxEs, start = - pi/6, end = pi/6), linetype = 1, show.legend = FALSE) +     
          #plot lines of equal lodes parameter and label them
          geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, alpha = lineAlpha), linetype = 1, data = lines, show.legend = FALSE) +
          geom_text(aes(x=x2, y=y2,label = label, vjust = -0.5, family = "Helvetica"), size = 2.75, data = lines) +
          #geom_text(aes(x=x, y=y, label = labels, angle = angle, family = "Times"), size = 2.75, data = axisLabels)+
          #plot data
          geom_point(aes(x=lodesBoots, y=octBoots), data = xyBoots, size = 1, alpha = 0.5)+
          geom_point(aes(x=lodes, y=oct), data = xy, size = 2, shape = 21, fill = "white")+
          geom_point(aes(x=meanLode, y=meanOct), data = xyMean, size = 4, color = "blue")+
          theme_void() +
          theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
     
     return(HsuPlot)
}


computeStationMeans <- function(amsData, stationList) {
     
     amsDataMeans <- data.frame()
     
     for (i in 1:length(stationList)){
          mean <- stationMean(stationList[[i]],amsData)
          amsDataMeans <- rbind(amsDataMeans, mean)
     }
     return(amsDataMeans)
}

PlotJellinekConfidenceIntervals <- function(stationName, amsData){
     data <- subset(amsData, grepl(paste(stationName), amsData$station))
     N = length(data$Name)
     dataMean <- ellMean(data$vector)
     #####     
     if (N > 5){
          
          
          vectorJ <- list()
          
          for (i in 1:length(data$Name)){
               vectorJ[[i]] <- c(data$tensor[[i]][1,1],data$tensor[[i]][2,2],data$tensor[[i]][3,3], data$tensor[[i]][1,2], data$tensor[[i]][2,3], data$tensor[[i]][1,3])
          }
          
          data$vectorJ <- vectorJ
          
          
          # Mean of ellipsoids
          
          
          
          vectorJ_mean <- list(c(dataMean$tensor[1,1], dataMean$tensor[2,2],dataMean$tensor[3,3],dataMean$tensor[1,2],dataMean$tensor[2,3],dataMean$tensor[1,3]))
          dataMean$vectorJ <- vectorJ_mean
          
          eigs <- eigen(dataMean$tensor)
          
          eigVecs <- eigs$vectors
          eigVals <- eigs$values
          
          
          # Compute Covariance matrix
          
          C <- matrix(nrow=6,ncol=6)
          
          for (j in 1:6) {
               for (k in 1:6){
                    tempSum = 0
                    for (l in 1:N){
                         tempSum = tempSum + (data$vectorJ[[l]][j] - dataMean$vectorJ[[1]][j]) * (data$vectorJ[[l]][k] - dataMean$vectorJ[[1]][k])
                    }
                    tempSum <- tempSum / N
                    
                    C[j,k] <- tempSum
               }
          }
          
          # Define the G Matrix
          Gcol1 <- c(eigVecs[1,1]^2, 
                     eigVecs[1,2]^2, 
                     eigVecs[1,3]^2, 
                     eigVecs[1,1]*eigVecs[1,2], 
                     eigVecs[1,2]*eigVecs[1,3], 
                     eigVecs[1,3]*eigVecs[1,1]) 
          Gcol2 <- c(eigVecs[2,1]^2, 
                     eigVecs[2,2]^2, 
                     eigVecs[2,3]^2, 
                     eigVecs[2,1]*eigVecs[2,2], 
                     eigVecs[2,2]*eigVecs[2,3], 
                     eigVecs[2,3]*eigVecs[2,1])
          Gcol3 <- c(eigVecs[3,1]^2,
                     eigVecs[3,2]^2, 
                     eigVecs[3,3]^2, 
                     eigVecs[3,1]*eigVecs[3,2], 
                     eigVecs[3,2]*eigVecs[3,3], 
                     eigVecs[3,3]*eigVecs[3,1]) 
          
          Gcol4 <- c(2*eigVecs[1,1]*eigVecs[2,1], 
                     2*eigVecs[1,2]*eigVecs[2,2], 
                     2*eigVecs[1,3]*eigVecs[2,3], 
                     eigVecs[1,1]*eigVecs[2,2] + eigVecs[2,1]*eigVecs[1,2],  
                     eigVecs[1,2]*eigVecs[2,3] + eigVecs[2,2]*eigVecs[1,3], 
                     eigVecs[1,3]*eigVecs[2,1] + eigVecs[2,3]*eigVecs[1,1])
          
          Gcol5 <- c(2*eigVecs[2,1]*eigVecs[3,1], 
                     2*eigVecs[2,2]*eigVecs[3,2], 
                     2*eigVecs[2,3]*eigVecs[3,3], 
                     eigVecs[2,1]*eigVecs[3,2] + eigVecs[3,1]*eigVecs[2,2],  
                     eigVecs[2,2]*eigVecs[3,3] + eigVecs[3,2]*eigVecs[2,3], 
                     eigVecs[2,3]*eigVecs[3,1] + eigVecs[3,3]*eigVecs[2,1])           
          
          Gcol6 <- c(2*eigVecs[3,1]*eigVecs[1,1], 
                     2*eigVecs[3,2]*eigVecs[1,2], 
                     2*eigVecs[3,3]*eigVecs[1,3], 
                     eigVecs[3,1]*eigVecs[1,2] + eigVecs[1,1]*eigVecs[3,2],  
                     eigVecs[3,2]*eigVecs[1,3] + eigVecs[1,2]*eigVecs[3,3], 
                     eigVecs[3,3]*eigVecs[1,1] + eigVecs[1,3]*eigVecs[3,1])    
          
          
          G <- cbind(Gcol1, Gcol2, Gcol3, Gcol4, Gcol5, Gcol6)
          
          # Rotate Covariance matrix into eigen vector coordinate system
          
          Cprime = G %*% C %*% t(G)
          
          
          
          # Find the 2D covarience matrix for each eigenvector
          i = 1
          j = 2
          k = 3
          
          W3 <- c(Cprime[i + 3, i + 3]/((eigVals[i]-eigVals[j])^2), 
                  Cprime[i + 3, k + 3]/((eigVals[i]-eigVals[j])*(eigVals[i]-eigVals[k])),
                  Cprime[i + 3, k + 3]/((eigVals[i]-eigVals[j])*(eigVals[i]-eigVals[k])),
                  Cprime[k + 3, k + 3]/((eigVals[i]-eigVals[k])^2))
          
          dim(W3) <- c(2,2)
          
          
          i = 3
          j = 1
          k = 2
          
          W1 <- c(Cprime[i + 3, i + 3]/((eigVals[i]-eigVals[j])^2), 
                  Cprime[i + 3, k + 3]/((eigVals[i]-eigVals[j])*(eigVals[i]-eigVals[k])),
                  Cprime[i + 3, k + 3]/((eigVals[i]-eigVals[j])*(eigVals[i]-eigVals[k])),
                  Cprime[k + 3, k + 3]/((eigVals[i]-eigVals[k])^2))
          
          dim(W1) <- c(2,2)
          
          
          i = 2
          j = 3
          k = 1
          
          W2 <- c(Cprime[i + 3, i + 3]/((eigVals[i]-eigVals[j])^2), 
                  Cprime[i + 3, k + 3]/((eigVals[i]-eigVals[j])*(eigVals[i]-eigVals[k])),
                  Cprime[i + 3, k + 3]/((eigVals[i]-eigVals[j])*(eigVals[i]-eigVals[k])),
                  Cprime[k + 3, k + 3]/((eigVals[i]-eigVals[k])^2))
          
          dim(W2) <- c(2,2)
          
          
          
          
          
          # Compute the degrees of the confidence ellipses
          
          W1eigs <- eigen(W1, symmetric = TRUE)
          W2eigs <- eigen(W2, symmetric = TRUE)
          W3eigs <- eigen(W3, symmetric = TRUE)
          
          W1eigVals <- W1eigs$values
          W2eigVals <- W2eigs$values
          W3eigVals <- W3eigs$values
          
          
          xy1 <- data.frame(t(W1eigs$vectors)) 
          xy2 <- data.frame(t(W2eigs$vectors))
          xy3 <- data.frame(t(W3eigs$vectors))
          
          
          Fstat <-(2*(N-1)/((N)*((N)-2)))*qf(.90,df1=2, df2=N-2) 
          
          JellEll_W1 <- c( ( atan( (Fstat*W1eigVals[1])^(1/2) ) ), ( atan( (Fstat*W1eigVals[2])^(1/2) ) ))
          rad2deg(JellEll_W1)
          JellEll_W2 <- c(( atan( (Fstat*W2eigVals[1])^(1/2) ) ), ( atan( (Fstat*W2eigVals[2])^(1/2) ) ))
          rad2deg(JellEll_W2)
          JellEll_W3 <- c(( atan( (Fstat*W3eigVals[1])^(1/2) ) ), ( atan( (Fstat*W3eigVals[2])^(1/2) ) ))
          rad2deg(JellEll_W3)
          
          
          
          
          # Plot the confidence regions for each. 
          numPoints = 1000
          
          #Confidence ellipse for 1st axis of mean 
          rot1 <- rbind((dataMean$rotation[3,]),(dataMean$rotation[1,]),(dataMean$rotation[2,]))
          result1 <- list()
          eig1 <- W1eigs
          semiaxes1 <- JellEll_W1
          circle <- lapply(0:numPoints, function(s) c(cos(s * 2 * pi / numPoints), sin(s * 2 * pi / numPoints)))
          vs1 <- lapply(circle, function(v) as.numeric((eig1$vectors) %*% (semiaxes1 * v)))
          result1$points <- lapply(vs1,rayPointFromTangentVector,rot1)
          
          
          #Confidence ellipse for 2nd axis of mean
          rot2 <- rbind((dataMean$rotation[2,]), (dataMean$rotation[3,]),(dataMean$rotation[1,]))
          
          result2 <- list()
          eig2 <- W2eigs
          semiaxes2 <- JellEll_W2
          circle <- lapply(0:numPoints, function(s) c(cos(s * 2 * pi / numPoints), sin(s * 2 * pi / numPoints)))
          vs2 <- lapply(circle, function(v) as.numeric((eig2$vectors) %*% (semiaxes2 * v)))
          result2$points <- lapply(vs2,rayPointFromTangentVector,rot2)
          
          #Confidence ellipse for 3rd axis of mean
          rot3 <- rbind((dataMean$rotation[1,]),(dataMean$rotation[2,]),(dataMean$rotation[3,]))
          result3 <- list()
          eig3 <- W3eigs
          semiaxes3 <- JellEll_W3
          circle <- lapply(0:numPoints, function(s) c(cos(s * 2 * pi / numPoints), sin(s * 2 * pi / numPoints)))
          vs3 <- lapply(circle, function(v) as.numeric((eig3$vectors) %*% (semiaxes3 * v)))
          result3$points <- lapply(vs3, rayPointFromTangentVector, rot3)
          
          
          
          K1 = data.frame(t(sapply(data$rotation, function(s) equalAreaProjection(lower(s[1,])))))
          K2 = data.frame(t(sapply(data$rotation, function(s) equalAreaProjection(lower(s[2,])))))
          K3 = data.frame(t(sapply(data$rotation, function(s) equalAreaProjection(lower(s[3,])))))
          Kmean1 = data.frame(t(equalAreaProjection(lower(dataMean$rotation[3,]))))
          Kmean2 = data.frame(t(equalAreaProjection(lower(dataMean$rotation[2,]))))
          Kmean3 = data.frame(t(equalAreaProjection(lower(dataMean$rotation[1,]))))
          
          
          
          points1 <- data.frame(t(sapply(result1$points, function(s) equalAreaProjection(lower(s)))))
          points2 <- data.frame(t(sapply(result2$points, function(s) equalAreaProjection(lower(s)))))
          points3 <- data.frame(t(sapply(result3$points, function(s) equalAreaProjection(lower(s)))))
          
          
          xy = data.frame(rbind(K1,K2,K3))
          xy = data.frame(xy, "axis"=rep(c(1,2,3), each=N))
          meanXy = data.frame(rbind(Kmean1, Kmean2, Kmean3))
          meanXy = data.frame(meanXy, "axis" = rep(c(1,2,3), each = 1))
          
          shapeVector = c(22,24,21)
          names(shapeVector) = c("1","2","3")
          
          colorVector = c("red","green","blue")
          names(colorVector) = c("1","2","3")
          
          
          ellBoots <- ellBootstrapInference(data$vector, numBoots = 1000)
          ells095 <- ellHighEllipsoidVectors(ellBoots$covarInv, ellBoots$center, (ellBoots$q095), numSamples=7)
          
          rotations095 <- list()
          
          for(i in 1:length(ells095)) {
               rotation <- ellEllipsoidFromVector(ells095[[i]])
               rotations095$rotation[[i]] = rotation$rotation
          }
          
          
          #PLOT: Bootstrapped ellipsoid vectors
          #NOTE: for real AMS data, should be S[3,], S[2,], S[1,]. For some reason it is flipped. 
          K1Boots = data.frame(t(sapply(rotations095$rotation, function(s) equalAreaProjection(lower(s[3,])))))
          K2Boots = data.frame(t(sapply(rotations095$rotation, function(s) equalAreaProjection(lower(s[2,])))))
          K3Boots = data.frame(t(sapply(rotations095$rotation, function(s) equalAreaProjection(lower(s[1,])))))
          xyBoots = data.frame(rbind(K1Boots,K2Boots,K3Boots))
          xyBoots = data.frame(xyBoots, "axis"=rep(c(1,2,3), each=nrow(K1Boots)))
          
          
          ggplot() +
               geom_point(data=xyBoots, aes(x=X1, y=X2, color = as.factor(axis),alpha = .005),show.legend = FALSE) +
               geom_point(data=points1, aes(x=X1, y=X2), color = "black", inherit.aes = FALSE, size=.01) +
               geom_point(data=points2, aes(x=X1, y=X2), color = "black", size = .01) +
               geom_point(data=points3, aes(x=X1, y=X2), color = "black", size= .01) +
               geom_point(data=xy, aes(x=X1, y=X2, shape = as.factor(axis), fill=as.factor(axis)), color = "black", size = 3, show.legend = FALSE) +
               geom_point(data=meanXy, aes(x=X1, y=X2, shape = as.factor(axis)), fill="black", color = "black", size = 5) +
               geom_circle(aes(x0=0,y0=0,r=sqrt(2)), inherit.aes = FALSE) +
               coord_fixed(ratio = 1) +
               scale_shape_manual(values = shapeVector) +
               scale_color_manual(values = colorVector) +
               theme_void()+
               ggtitle(paste(data$Name[[1]]))
          
          
     } 
     else {
          
          K1 = data.frame(t(sapply(data$rotation, function(s) equalAreaProjection(lower(s[1,])))))
          K2 = data.frame(t(sapply(data$rotation, function(s) equalAreaProjection(lower(s[2,])))))
          K3 = data.frame(t(sapply(data$rotation, function(s) equalAreaProjection(lower(s[3,])))))
          Kmean1 = data.frame(t(equalAreaProjection(lower(dataMean$rotation[3,]))))
          Kmean2 = data.frame(t(equalAreaProjection(lower(dataMean$rotation[2,]))))
          Kmean3 = data.frame(t(equalAreaProjection(lower(dataMean$rotation[1,]))))
          
          
          
          xy = data.frame(rbind(K1,K2,K3))
          xy = data.frame(xy, "axis"=rep(c(1,2,3), each=N))
          meanXy = data.frame(rbind(Kmean1, Kmean2, Kmean3))
          meanXy = data.frame(meanXy, "axis" = rep(c(1,2,3), each = 1))
          
          shapeVector = c(22,24,21)
          names(shapeVector) = c("1","2","3")
          
          colorVector = c("red","green","blue")
          names(colorVector) = c("1","2","3")
          
     
          
          
          ggplot() +
               geom_point(data=xy, aes(x=X1, y=X2, shape = as.factor(axis), fill=as.factor(axis)), color = "black", size = 3, show.legend = FALSE) +
               geom_point(data=meanXy, aes(x=X1, y=X2, shape = as.factor(axis)), fill="black", color = "black", size = 5) +
               geom_circle(aes(x0=0,y0=0,r=sqrt(2)), inherit.aes = FALSE) +
               coord_fixed(ratio = 1) +
               scale_shape_manual(values = shapeVector) +
               scale_color_manual(values = colorVector) +
               theme_void()+
               ggtitle(paste(data$Name[[1]]))
          
     }
}

plotAllJellAndBootstraps <- function(amsData, stationList) {
     bootstrappedData <- data.frame()
     bootPlots <- list()
     notBootstrapped <- list()
     count = 1
     count2 = 1
     for (i in 1:length(stationList)){
          print(paste(i, " of ", length(stationList)))
          bootPlots[[i]] <- PlotJellinekConfidenceIntervals(stationList[[i]], amsData)
          
          
          
          # boots95 <- computeAMSBootStrap95(stationList[[i]], amsData)
          # if (typeof(boots95) != "character") {
          #      bootPlots[[count]] <- plotNetBoots(boots95)
          #      bootstrappedData <- rbind(bootstrappedData, boots95)
          #      count = count + 1
          # }
          # else {
          #      
          #      notBootstrapped[[count2]] <- stationList[[i]]
          #      count2 = count2+1
          # }
     }
     return(bootPlots)
}

plotShapeVersusAnisotropy <- function(amsData, gradient = "lodes") {
     lodes = sapply(amsData$logA, function(s) ellLodeNu(s))
     oct = sapply(amsData$logA, function(s) ellOctahedralShearStrain(s))
     Km = sapply(amsData$Km, function(s) s)
     lodesOct = data.frame(cbind(oct, lodes, Km))
     if (gradient == "lodes") {
      plot <- ggplot(data=lodesOct, aes(x=oct, y=lodes, color = lodes)) +
                    geom_point(show.legend = TRUE) +
                    theme_classic()+
                    coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
                    geom_abline(slope = 0, intercept = 0) +
                    xlab("anisotropy (Oct. Shear)") +
                    ylab("shape (Lodes)")+
                    ggtitle(label = "Shape and Anisotropy parameters", subtitle = "colored by Lodes shape parameter") +
                    scale_colour_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", guide = "colourbar")
      
     } else if (gradient == "oct") {
         plot <-  ggplot(data=lodesOct, aes(x=oct, y=lodes, color = oct)) +
               geom_point(show.legend = TRUE) +
               theme_classic()+
               coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               geom_abline(slope = 0, intercept = 0) +
               xlab("anisotropy (Oct. Shear)") +
               ylab("shape (Lodes)")+
               ggtitle(label = "Shape and Anisotropy parameters", subtitle = "colored by anisotropy degree (oct)") +
               scale_colour_gradient(low = "#ffffbf", high = "#4575b4", guide = "colourbar")
         
     } else if (gradient == "Km") {
          plot <- ggplot(data=lodesOct, aes(x=oct, y=lodes, color = log10(Km))) +
               geom_point(show.legend = TRUE) +
               theme_classic()+
               coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               geom_abline(slope = 0, intercept = 0) +
               xlab("anisotropy (Oct. Shear)") +
               ylab("shape (Lodes)")+
               ggtitle(label = "Shape and Anisotropy parameters", subtitle = "colored by bulk susceptibility (Km)") +
               scale_colour_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", midpoint = -3.2, guide = "colourbar")
     }
     
     else {
          plot <- ggplot(data=lodesOct, aes(x=oct, y=lodes)) +
               geom_point() +
               theme_classic()+
               coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               geom_abline(slope = 0, intercept = 0) +
               xlab("anisotropy (Oct. Shear)") +
               ylab("shape (Lodes)") +
               ggtitle(label = "Shape and Anisotropy parameters") 
               
     }
     
     plot
}

plotHistogramLodes <- function(amsData, binWidth) {
     lodes = (sapply(amsData$logA, function(s) ellLodeNu(s)))
     lodes = data.frame(lodes)
     ggplot() +
          geom_histogram(data = lodes, aes(x = lodes), color = "white", fill = "black", binwidth = binWidth) +
          theme_classic()
     }

plotShapeVersusKm <- function(amsData, gradient = "Km") {
     lodes = sapply(amsData$logA, function(s) ellLodeNu(s))
     oct = sapply(amsData$logA, function(s) ellOctahedralShearStrain(s))
     Km = sapply(amsData$Km, function(s) log10(s))
     lodesOct = data.frame(cbind(oct, lodes, Km))
     if (gradient == "lodes") {
          plot <- ggplot(data=lodesOct, aes(x=Km, y=lodes, color = lodes)) +
               geom_point(show.legend = TRUE) +
               theme_classic()+
               #coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               geom_abline(slope = 0, intercept = 0, color = "grey") +
               xlab("log10(Km)") +
               ylab("shape (Lodes)")+
               ggtitle(label = "Lodes Shape Parameter and Susceptibility", subtitle = "colored by Lodes shape parameter") +
               scale_colour_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", guide = "colourbar")
          
     } else if (gradient == "oct") {
          plot <-  ggplot(data=lodesOct, aes(x=Km, y=lodes, color = oct)) +
               geom_point(show.legend = TRUE) +
               theme_classic()+
               #coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               geom_abline(slope = 0, intercept = 0) +
               xlab("log10(Km)") +
               ylab("shape (Lodes)")+
               ggtitle(label = "Lodes Shape Parameter and Susceptibility", subtitle = "colored by anisotropy degree (oct)") +
               scale_colour_gradient(low = "#ffffbf", high = "#4575b4", guide = "colourbar")
          
     } else if (gradient == "Km") {
          plot <- ggplot(data=lodesOct, aes(x=Km, y=lodes, color = Km)) +
               geom_point(show.legend = TRUE) +
               theme_classic()+
               #coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               geom_abline(slope = 0, intercept = 0) +
               xlab("log10(Km)") +
               ylab("shape (Lodes)")+
               ggtitle(label = "Lodes Shape Parameter and Susceptibility", subtitle = "colored by bulk susceptibility (Km)") +
               scale_colour_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", midpoint = -3.2, guide = "colourbar")
     }
     
     else {
          plot <- ggplot(data=lodesOct, aes(x=Km, y=lodes)) +
               geom_point() +
               theme_classic()+
               #coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               geom_abline(slope = 0, intercept = 0) +
               ggtitle(label = "Lodes Shape Parameter and Susceptibility") +
               xlab("log10(Km)") +
               ylab("shape (Lodes)")
     }
     
     plot
}

plotOctVersusKm <- function(amsData, gradient = "Km") {
     lodes = sapply(amsData$logA, function(s) ellLodeNu(s))
     oct = sapply(amsData$logA, function(s) ellOctahedralShearStrain(s))
     Km = sapply(amsData$Km, function(s) log10(s))
     lodesOct = data.frame(cbind(oct, lodes, Km))
     if (gradient == "lodes") {
          plot <- ggplot(data=lodesOct, aes(x=Km, y=oct, color = lodes)) +
               geom_point() +
               theme_classic()+
               #coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               #geom_abline(slope = 0, intercept = 0, color = "grey") +
               xlab("log10(Km)") +
               ylab("anisotropy degree (oct)")+
               ggtitle(label = "Anisotropy degree and susceptibility", subtitle = "colored by Lodes shape parameter") +
               scale_colour_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", guide = "colourbar")
          
          
     } else if (gradient == "oct") {
        plot <- ggplot(data=lodesOct, aes(x=Km, y=oct, color = oct)) +
               geom_point() +
               theme_classic()+
               #coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               #geom_abline(slope = 0, intercept = 0) +
               xlab("log10(Km)") +
               ylab("anisotropy degree (oct)")+
               ggtitle(label = "Anisotropy degree and susceptibility", subtitle = "colored by anisotropy degree (oct)") +
               scale_colour_gradient(low = "#ffffbf", high = "#4575b4", guide = "colourbar")
          
     } else if (gradient == "Km") {
          plot <- ggplot(data=lodesOct, aes(x=Km, y=oct, color = Km)) +
               geom_point() +
               theme_classic()+
               #coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               #geom_abline(slope = 0, intercept = 0) +
               xlab("log10(Km)") +
               ylab("anisotropy degree (oct)")+
               ggtitle(label = "Anisotropy degree and susceptibility", subtitle = "colored by bulk susceptibility (Km)") +
               scale_colour_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", midpoint = -3.2, guide = "colourbar")
     }
     
     else {
          plot <- ggplot(data=lodesOct, aes(x=Km, y=oct)) +
               geom_point() +
               theme_classic()+
               #coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               #geom_abline(slope = 0, intercept = 0) +
               xlab("log10(Km)") +
               ylab("anisotropy degree (oct)") +
               ggtitle(label = "Anisotropy degree and susceptibility") 
     }
     
     plot
}



plotLodesVersusKm <- function(amsData, gradient = "Km") {
     lodes = sapply(amsData$logA, function(s) ellLodeNu(s))
     oct = sapply(amsData$logA, function(s) ellOctahedralShearStrain(s))
     Km = sapply(amsData$Km, function(s) log10(s))
     lodesOct = data.frame(cbind(oct, lodes, Km))
     if (gradient == "lodes") {
          plot <- ggplot(data=lodesOct, aes(x=Km, y=lodes, color = lodes)) +
               geom_point() +
               theme_classic()+
               #coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               #geom_abline(slope = 0, intercept = 0, color = "grey") +
               xlab("log10(Km)") +
               ylab("anisotropy degree (oct)")+
               ggtitle(label = "Anisotropy degree and susceptibility", subtitle = "colored by Lodes shape parameter") +
               scale_colour_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", guide = "colourbar")
          
          
     } else if (gradient == "oct") {
          plot <- ggplot(data=lodesOct, aes(x=Km, y=lodes, color = oct)) +
               geom_point() +
               theme_classic()+
               #coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               #geom_abline(slope = 0, intercept = 0) +
               xlab("log10(Km)") +
               ylab("anisotropy degree (oct)")+
               ggtitle(label = "Anisotropy degree and susceptibility", subtitle = "colored by anisotropy degree (oct)") +
               scale_colour_gradient(low = "#ffffbf", high = "#4575b4", guide = "colourbar")
          
     } else if (gradient == "Km") {
          plot <- ggplot(data=lodesOct, aes(x=Km, y=lodes, color = Km)) +
               geom_point() +
               theme_classic()+
               #coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               #geom_abline(slope = 0, intercept = 0) +
               xlab("log10(Km)") +
               ylab("anisotropy degree (oct)")+
               ggtitle(label = "Anisotropy degree and susceptibility", subtitle = "colored by bulk susceptibility (Km)") +
               scale_colour_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", midpoint = -3.2, guide = "colourbar")
     }
     
     else {
          plot <- ggplot(data=lodesOct, aes(x=Km, y=lodes)) +
               stat_density2d(aes(fill = ..level..),geom = "polygon", n = 100) +
               scale_fill_viridis() +
               geom_point(size = 1, alpha = 0.5) +
               theme_classic()+
               #coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
               #geom_abline(slope = 0, intercept = 0) +
               xlab("log10(Km)") +
               ylab("anisotropy degree (oct)") +
               ggtitle(label = "Anisotropy degree and susceptibility") +
               geom_smooth(method='lm',formula=y~x)
          }
     
     plot
}












plotMapByProperty <- function(amsData, property = "Km", geoMap) {
     
     
     if (property == "Km") {
          plot <- ggplot(geoMap["SUPERSUITE"]) +
               geom_sf(aes(fill = SUPERSUITE)) +
               coord_sf(ylim=c(-20.85,-21.5),xlim=c(119.75,120.4), expand = FALSE) +
               scale_fill_manual(values=c("Split Rock Supersuite" = "#F77D84", "Cleland Supersuite" = "#EEC0DA", "Emu Pool Supersuite" = "#E49DC5", "Tambina Supersuite" = "#DA78B1", "Callina Supersuite" = "#D557AC")) +
               geom_point(data=amsData, aes(x=easting, y=northing, color=log10(Km)), stroke = 0, size = 6) +
               scale_color_viridis() +
               #theme_bw() +
               theme_void() +
               ggtitle(label = "AMS station locations", subtitle = "colored by bulk susceptibility (Km)")
     }
     else if (property == "oct") {
          oct = sapply(amsData$logA, function(s) ellOctahedralShearStrain(s))
          amsData$oct <- oct
          
          plot <-  ggplot(geoMap["SUPERSUITE"]) +
               geom_sf(aes(fill = SUPERSUITE)) +
               coord_sf(ylim=c(-20.85,-21.5),xlim=c(119.75,120.4), expand = FALSE) +
               scale_fill_manual(values=c("Split Rock Supersuite" = "#F77D84", "Cleland Supersuite" = "#EEC0DA", "Emu Pool Supersuite" = "#E49DC5", "Tambina Supersuite" = "#DA78B1", "Callina Supersuite" = "#D557AC")) +
               
               geom_point(data=amsData, aes(x=easting, y=northing, color=oct), stroke = 0, size = 4) +
               scale_color_viridis(limits = c(0,0.2)) +
               theme_bw() +
               ggtitle(label = "AMS station locations", subtitle = "colored by anisotropy degree (oct)")
          
     }
     
     else if (property == "lodes") {
          lodes = sapply(amsData$logA, function(s) ellLodeNu(s))
          amsData$lodes <- lodes
          
          # plot <-  ggplot(geoMap["SUPERSUITE"]) +
          #      geom_sf(aes(fill = SUPERSUITE)) +
          #      coord_sf(ylim=c(-20.85,-21.5),xlim=c(119.75,120.4), expand = FALSE) +
          #      scale_fill_manual(values=c("Split Rock Supersuite" = "#F77D84", "Cleland Supersuite" = "#EEC0DA", "Emu Pool Supersuite" = "#E49DC5", "Tambina Supersuite" = "#DA78B1", "Callina Supersuite" = "#D557AC")) +
          #      geom_point(data=amsData, aes(x=easting, y=northing, color=lodes), stroke = 0, size = 4) +
          #      scale_color_viridis(limits = c(-1,1) ) +
          #      theme_bw() +
          #      ggtitle(label = "AMS station locations", subtitle = "colored by lodes shape parameter")
          # 
          plot <-  ggplot(geoMap["SUPERSUITE"]) +
               geom_sf(aes(fill = SUPERSUITE)) +
               coord_sf(ylim=c(-20.85,-21.5),xlim=c(119.75,120.4), expand = FALSE) +
               scale_fill_manual(values=c("Split Rock Supersuite" = "#F77D84", "Cleland Supersuite" = "#EEC0DA", "Emu Pool Supersuite" = "#E49DC5", "Tambina Supersuite" = "#DA78B1", "Callina Supersuite" = "#D557AC")) +
               geom_point(data=amsData, aes(x=easting, y=northing, color=lodes), stroke = 0, size = 6) +
               scale_color_viridis(limits = c(-1,1) ) +
               theme_void() +
               ggtitle(label = "AMS station locations", subtitle = "colored by lodes shape parameter")
          
     }
     
     else {
          print(paste0("'", property, "' is not a known property. Your choices are 'Km', 'oct', or 'lodes'"))
                
                
          plot <- ggplot(geoMap["SUPERSUITE"]) +
               geom_sf(aes(fill = SUPERSUITE)) +
               coord_sf(ylim=c(-20.85,-21.5),xlim=c(119.75,120.4), expand = FALSE) +
               scale_fill_manual(values=c("Split Rock Supersuite" = "#F77D84", "Cleland Supersuite" = "#EEC0DA", "Emu Pool Supersuite" = "#E49DC5", "Tambina Supersuite" = "#DA78B1", "Callina Supersuite" = "#D557AC")) +
               
               geom_point(data=amsData, aes(x=easting, y=northing), stroke = 0, size = 4, show.legend = FALSE) +
                     theme_bw() +
                     ggtitle(label = "AMS station locations")
          
                
     }
     
     plot

}


plotMapPlots <- function(dataSet, geoMap, path) {
     
     imgList <- list()
     
     for (i in 1:length(dataSet$Name)) {
          img <- load.image(paste0(path, "/",dataSet$Name[[i]],".png"))
          imgList[i] <- rasterGrob(img, interpolate=FALSE)
     }
     
     
     dataSet$plot  <- imgList
     
     plot <- ggplot(geoMap["SUPERSUITE"]) +
          geom_sf(aes(fill = SUPERSUITE)) +
          coord_sf(ylim=c(-20.85,-21.5),xlim=c(119.75,120.4), expand = FALSE) +
          scale_fill_manual(values=c("Split Rock Supersuite" = "#F77D84", "Cleland Supersuite" = "#EEC0DA", "Emu Pool Supersuite" = "#E49DC5", "Tambina Supersuite" = "#DA78B1", "Callina Supersuite" = "#D557AC")) +
          #geom_point(data = dataSet, aes(x = as.numeric(dataSet$easting), y = as.numeric(dataSet$northing), color = dip), size = 5) + 
          #scale_color_viridis(limits = c(0,90) ) +
          geom_blank() +
          mapply(function(xx, yy, imgs) {
               annotation_raster(imgs, xmin=xx-0.005, xmax=xx+.005, ymin=yy-0.0025, ymax=yy+0.0025)},
               as.numeric(dataSet$easting), as.numeric(dataSet$northing), (dataSet$plot)) +
          theme_bw() +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank())
     
     plot
}


plotMapStrikesDips <- function(dataSet, geoMap, onlyStat = FALSE, xlims = c(119.75,120.4), ylims = c(-20.85,-21.5) ){
     
     #read file
     img <- load.image("StrikeDipFol.png")
     imgSad <- load.image("sadFace.png")
     
     sd <- lapply(dataSet$rotation, function(s) geoStrikeDipDegFromCartesian(s[3,]))
     
     imgList <- list()
     
     for (i in 1:length(dataSet$Name)) {
          
          if (onlyStat == FALSE) {
               rotImage <- imrotate(img, sd[[i]][1])
               rotImage2 <- resize(rotImage, size_x = 130L, size_y= 130L, interpolation_type = 0, centering_x = 0.5, centering_y = 0.5)
               imgList[i] <- rasterGrob(rotImage2, interpolate=FALSE)
               imgList[[i]][imgList[[i]] != "#1D0047"] <- "#FFFFFF00"
               imgList[[i]][imgList[[i]] == "#1D0047"] <- "#000000FF"
          }
          
          if (onlyStat == TRUE) {
               if (dataSet$specimens[[i]] < 5) {
                    rotImage <- imrotate(imgSad, sd[[i]][1])
                    rotImage2 <- resize(rotImage, size_x = 130L, size_y= 130L, interpolation_type = 0, centering_x = 0.5, centering_y = 0.5)
                    imgList[i] <- rasterGrob(rotImage2, interpolate=FALSE)
                    imgList[[i]][imgList[[i]] != "#1D0047"] <- "#FFFFFF00"
                    imgList[[i]][imgList[[i]] == "#1D0047"] <- "#000000FF"
               }
               else {
                    rotImage <- imrotate(img, sd[[i]][1])
                    rotImage2 <- resize(rotImage, size_x = 130L, size_y= 130L, interpolation_type = 0, centering_x = 0.5, centering_y = 0.5)
                    imgList[i] <- rasterGrob(rotImage2, interpolate=FALSE)
                    imgList[[i]][imgList[[i]] != "#1D0047"] <- "#FFFFFF00"
                    imgList[[i]][imgList[[i]] == "#1D0047"] <- "#000000FF"   
               }
          }
     }
     
     #dataSet$sdSymbol <- replicate(length(dataSet$Name), 1)
     
     dataSet$sdSymbol  <- lapply(imgList, function(s) s)
     dataSet$strike <- sapply(sd, function(s) s[1])
     dataSet$dip <- sapply(sd, function(s) s[2])
     
     plot <- ggplot(geoMap["SUPERSUITE"]) +
          geom_sf(aes(fill = SUPERSUITE)) +
          coord_sf(ylim=ylims,xlim=xlims, expand = FALSE) +
          scale_fill_manual(values=c("Split Rock Supersuite" = "#F77D84", "Cleland Supersuite" = "#EEC0DA", "Emu Pool Supersuite" = "#E49DC5", "Tambina Supersuite" = "#DA78B1", "Callina Supersuite" = "#D557AC")) +
          geom_point(data = dataSet, aes(x = as.numeric(dataSet$easting), y = as.numeric(dataSet$northing), color = dip), size = 6) + 
          scale_color_viridis(limits = c(0,90) ) +
          geom_blank() +
          mapply(function(xx, yy, imgs) {
               annotation_raster(imgs, xmin=xx-.007, xmax=xx+.007, ymin=yy-0.007, ymax=yy+0.007)},
               as.numeric(dataSet$easting), as.numeric(dataSet$northing), (dataSet$sdSymbol)) +
          #theme_bw() +
          theme_void() +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank())
     
     plot
}
     
#      
#      #PLOT: Lodes Parameter as a function of Octahedral shear strain
#      lodes = sapply(data$logA, function(s) ellLodeNu(s))
#      oct = sapply(data$logA, function(s) ellOctahedralShearStrain(s))
#      lodesOct = data.frame(cbind(oct, lodes))
#      
#      ggplot(data=lodesOct, aes(x=oct, y=lodes)) +
#           geom_point() +
#           theme_classic()+
#           coord_cartesian(ylim = c(-1,1), xlim = c(0,.3)) +
#           geom_abline(slope = 0, intercept = 0) +
#           xlab("anisotropy (Oct. Shear)") +
#           ylab("shape (Lodes)")
#      
# }
# allPlots[[1]]
# }

plotMapTrendsPlunges <- function(dataSet, geoMap, xlims = c(119.75,120.4), ylims = c(-20.85,-21.5), size = "normal" ){
     
     #read file
     Lin0to10 <- load.image("Symbol_Lin0_10.png")
     Lin10to20 <- load.image("Symbol_Lin10_20.png")
     Lin20to30 <- load.image("Symbol_Lin20_30.png")
     Lin30to40 <- load.image("Symbol_Lin30_40.png")
     Lin40to50 <- load.image("Symbol_Lin40_50.png")
     Lin50to60 <- load.image("Symbol_Lin50_60.png")
     Lin60to70 <- load.image("Symbol_Lin60_70.png")
     Lin70to80 <- load.image("Symbol_Lin70_80.png")
     Lin80to90 <- load.image("Symbol_Lin80_90.png")
     
     
     tp <- lapply(dataSet$rotation, function(s) geoTrendPlungeDegFromCartesian(lower(s[1,])))
     
     imgList <- list()
     
     for (i in 1:length(dataSet$Name)) {
          if (tp[[i]][2] <= 10) {
               rotImage <- imrotate(Lin0to10, tp[[i]][1])
               rotImage2 <- resize(rotImage, size_x = 400L, size_y= 400L, interpolation_type = 0, centering_x = 0.5, centering_y = 0.5)
               imgList[i] <- rasterGrob(rotImage2, interpolate=FALSE)
               imgList[[i]][imgList[[i]] != "#1D0047"] <- "#FFFFFF00"
               imgList[[i]][imgList[[i]] == "#1D0047"] <- "#000000FF"  
          }
          
          else if (tp[[i]][2] > 10 && tp[[i]][2] <= 20) {
               rotImage <- imrotate(Lin10to20, tp[[i]][1])
               rotImage2 <- resize(rotImage, size_x = 400L, size_y= 400L, interpolation_type = 0, centering_x = 0.5, centering_y = 0.5)
               imgList[i] <- rasterGrob(rotImage2, interpolate=FALSE)
               imgList[[i]][imgList[[i]] != "#1D0047"] <- "#FFFFFF00"
               imgList[[i]][imgList[[i]] == "#1D0047"] <- "#000000FF"  
          }
          
          else if (tp[[i]][2] > 20 && tp[[i]][2] <= 30) {
               rotImage <- imrotate(Lin20to30, tp[[i]][1])
               rotImage2 <- resize(rotImage, size_x = 400L, size_y= 400L, interpolation_type = 0, centering_x = 0.5, centering_y = 0.5)
               imgList[i] <- rasterGrob(rotImage2, interpolate=FALSE)
               imgList[[i]][imgList[[i]] != "#1D0047"] <- "#FFFFFF00"
               imgList[[i]][imgList[[i]] == "#1D0047"] <- "#000000FF"  
          }
          
          else if (tp[[i]][2] > 30 && tp[[i]][2] <= 40) {
               rotImage <- imrotate(Lin30to40, tp[[i]][1])
               rotImage2 <- resize(rotImage, size_x = 400L, size_y= 400L, interpolation_type = 0, centering_x = 0.5, centering_y = 0.5)
               imgList[i] <- rasterGrob(rotImage2, interpolate=FALSE)
               imgList[[i]][imgList[[i]] != "#1D0047"] <- "#FFFFFF00"
               imgList[[i]][imgList[[i]] == "#1D0047"] <- "#000000FF"  
          }
          
          else if (tp[[i]][2] > 40 && tp[[i]][2] <= 50) {
               rotImage <- imrotate(Lin40to50, tp[[i]][1])
               rotImage2 <- resize(rotImage, size_x = 400L, size_y= 400L, interpolation_type = 0, centering_x = 0.5, centering_y = 0.5)
               imgList[i] <- rasterGrob(rotImage2, interpolate=FALSE)
               imgList[[i]][imgList[[i]] != "#1D0047"] <- "#FFFFFF00"
               imgList[[i]][imgList[[i]] == "#1D0047"] <- "#000000FF"  
          }
          
          else if (tp[[i]][2] > 50 && tp[[i]][2] <= 60) {
               rotImage <- imrotate(Lin50to60, tp[[i]][1])
               rotImage2 <- resize(rotImage, size_x = 400L, size_y= 400L, interpolation_type = 0, centering_x = 0.5, centering_y = 0.5)
               imgList[i] <- rasterGrob(rotImage2, interpolate=FALSE)
               imgList[[i]][imgList[[i]] != "#1D0047"] <- "#FFFFFF00"
               imgList[[i]][imgList[[i]] == "#1D0047"] <- "#000000FF"  
          }
          
          else if (tp[[i]][2] > 60 && tp[[i]][2] <= 70) {
               rotImage <- imrotate(Lin60to70, tp[[i]][1])
               rotImage2 <- resize(rotImage, size_x = 400L, size_y= 400L, interpolation_type = 0, centering_x = 0.5, centering_y = 0.5)
               imgList[i] <- rasterGrob(rotImage2, interpolate=FALSE)
               imgList[[i]][imgList[[i]] != "#1D0047"] <- "#FFFFFF00"
               imgList[[i]][imgList[[i]] == "#1D0047"] <- "#000000FF"  
          }
          
          else if (tp[[i]][2] > 70 && tp[[i]][2] <= 80) {
               rotImage <- imrotate(Lin70to80, tp[[i]][1])
               rotImage2 <- resize(rotImage, size_x = 400L, size_y= 400L, interpolation_type = 0, centering_x = 0.5, centering_y = 0.5)
               imgList[i] <- rasterGrob(rotImage2, interpolate=FALSE)
               imgList[[i]][imgList[[i]] != "#1D0047"] <- "#FFFFFF00"
               imgList[[i]][imgList[[i]] == "#1D0047"] <- "#000000FF"  
          }
          
          else if (tp[[i]][2] > 80 && tp[[i]][2] <= 90) {
               rotImage <- imrotate(Lin80to90, tp[[i]][1])
               rotImage2 <- resize(rotImage, size_x = 400L, size_y= 400L, interpolation_type = 0, centering_x = 0.5, centering_y = 0.5)
               imgList[i] <- rasterGrob(rotImage2, interpolate=FALSE)
               imgList[[i]][imgList[[i]] != "#1D0047"] <- "#FFFFFF00"
               imgList[[i]][imgList[[i]] == "#1D0047"] <- "#000000FF"  
          }
          else {}
     }
     
     dataSet$tpSymbol  <- imgList
     dataSet$trend <- sapply(tp, function(s) s[1])
     dataSet$plunge <- sapply(tp, function(s) s[2])
     
     if (size == "small") {
     plot <- ggplot(geoMap["SUPERSUITE"]) +
          geom_sf(aes(fill = SUPERSUITE)) +
          coord_sf(ylim=ylims,xlim=xlims, expand = FALSE) +
          scale_fill_manual(values=c("Split Rock Supersuite" = "#F77D84", "Cleland Supersuite" = "#EEC0DA", "Emu Pool Supersuite" = "#E49DC5", "Tambina Supersuite" = "#DA78B1", "Callina Supersuite" = "#D557AC")) +
          geom_point(data = dataSet, aes(x = as.numeric(dataSet$easting), y = as.numeric(dataSet$northing), color = plunge), size = 6) + 
          scale_color_viridis(limits = c(0,90) ) +
          geom_blank() +
          mapply(function(xx, yy, imgs) {
               annotation_raster(imgs, xmin=xx-.02, xmax=xx+.02, ymin=yy-0.02, ymax=yy+0.02)},
               as.numeric(dataSet$easting), as.numeric(dataSet$northing), (dataSet$tpSymbol)) +
          theme_void() +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank())
     } else {
          plot <- ggplot(geoMap["SUPERSUITE"]) +
               geom_sf(aes(fill = SUPERSUITE)) +
               coord_sf(ylim=ylims,xlim=xlims, expand = FALSE) +
               scale_fill_manual(values=c("Split Rock Supersuite" = "#F77D84", "Cleland Supersuite" = "#EEC0DA", "Emu Pool Supersuite" = "#E49DC5", "Tambina Supersuite" = "#DA78B1", "Callina Supersuite" = "#D557AC")) +
               geom_point(data = dataSet, aes(x = as.numeric(dataSet$easting), y = as.numeric(dataSet$northing), color = plunge), size = 6) + 
               scale_color_viridis(limits = c(0,90) ) +
               geom_blank() +
               mapply(function(xx, yy, imgs) {
                    annotation_raster(imgs, xmin=xx-.05, xmax=xx+.05, ymin=yy-0.05, ymax=yy+0.05)},
                    as.numeric(dataSet$easting), as.numeric(dataSet$northing), (dataSet$tpSymbol)) +
               theme_void() +
               theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank())
     }
     
     plot
}



#Jelinek's "Corrected" anisotropy factor, P'
calcJelP <- function(logA){
     meanLogA <- (logA[1] + logA[2] + logA[3])/3
     JelP <- exp(
          (2*((logA[1] - meanLogA)^2 +(logA[2] - meanLogA)^2 + (logA[2] - meanLogA)^2))^(1/2)
     )
     return(JelP)
}

#Jelinek's shape factor, T
calcJelT <- function(logA) {
     JelT <- (2 * logA[2] - logA[1] - logA[3])/(logA[1]-logA[3])
     return(JelT)
}

plotShapeVersusAnisotropyBootstraps <- function(stationName, amsData) {
     amsStationData <- subset(amsData, grepl(paste(stationName), amsData$station))
     
     bootsFrame <- computeAMSBootStrap95(stationName, amsData) 
     meanData <- stationMean(stationName, amsData)
     # logAs <- list()
     # for(i in 1:length(ellBoots$bootstraps)) {
     #                ellEllipsoids <- ellEllipsoidFromVector(ellBoots$bootstraps[[i]])
     #                logAtemp <- c(ellEllipsoids$logA[3], ellEllipsoids$logA[2], ellEllipsoids$logA[1])
     #                logAs[[i]] <- logAtemp
     #           }
     
     JelT <- sapply(amsStationData$logA, function(s) calcJelT(s))
     JelP <- sapply(amsStationData$logA, function(s) ellJelinekP(s))
     JelTboots <- sapply(bootsFrame$logA, function(s) calcJelT(s))
     JelPboots <- sapply(bootsFrame$logA, function(s) ellJelinekP(s))
     JelTmean <- sapply(meanData$logA, function(s) calcJelT(s))
     JelPmean <- sapply(meanData$logA, function(s) ellJelinekP(s))
     
     PandT <- data.frame(cbind(JelP, JelT))
     names(PandT) <- c("jelP","jelT")
     
     PandTboots <- data.frame(cbind(JelPboots,JelTboots))
     names(PandTboots) = c("jelP","jelT")
     
     PandTmean <- data.frame(cbind(JelPmean,JelTmean))
     names(PandTmean) = c("jelP","jelT")
     
     ggplot() +
          geom_point(data=PandTboots, aes(x=jelP, y=jelT), alpha = .05) +
          geom_point(data=PandT, aes(x=jelP, y=jelT), shape = 21, fill = "white", color = "black", size = 3) +
          geom_point(data=PandTmean, aes(x=jelP, y=jelT), shape = 23, fill = "yellow", color = "black", size = 5) +
          theme_classic() +
          coord_cartesian(ylim = c(-1,1), xlim = c(1,1.3)) +
          geom_abline(slope = 0, intercept = 0) +
          xlab("Jelinek P' anisotropy factor") +
          ylab("Jelinek T shape factor") +
          ggtitle(paste("Shape and Anisotropy of", stationName))
     
     
}










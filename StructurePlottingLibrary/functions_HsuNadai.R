# These functions are dependent on GGPLOT2 and GGFORCE packages. make sure both are installed and loaded!


####This is a function to plot a blank Hsu-Nadai plot with various options.
#' @param plotArcs boolean. If TRUE, curved arcs of equal octaheral shear strain will be plotted. 
#' @param numArcs number of arcs to plot. An integer. 
#' @param plotLines boolean. If TRUE, light grey lines of equal lodes parameter will be plotted.  
#' @param maxOct the max octahedral shear strain to plot. a positive number.
plotBlankHsuNadaiPlot = function(plotArcs = TRUE, numArcs = 4, plotLines = TRUE, maxOct = 3) {
     # Make data frame of Hsu plot arcs of equal oct. shear strain
     # Change arcInterval to have more or fewer arcs display on the Hsu plot
     arcInterval = maxOct / numArcs
     
     arcs <- data.frame(
          start = rep(-pi/6, maxOct/arcInterval + 1),
          end = rep(pi/6, maxOct/arcInterval + 1),
          r = seq(from = 0, to = maxOct , by = arcInterval),
          textX = seq(from = 0*cos(pi/3), to = maxOct*cos(pi/3), by=arcInterval * cos(pi/3)),
          textY = seq(from = 0*sin(pi/3), to = maxOct*sin(pi/3), by=arcInterval * sin(pi/3))
     )
     
     lines <- data.frame(x1 = c(0,0,0), 
                         x2 = c(maxOct*cos(pi/3),maxOct*cos(pi/2), maxOct*cos(2/3 * pi)),
                         y1 = c(0,0,0), 
                         y2 = c(maxOct*sin(pi/3),maxOct*sin(pi/2), maxOct*sin(2/3 * pi)),
                         label = c(1,0, -1)
     )
     
     subLines <- data.frame(x1 = c(0,0,0,0,0,0), 
                            x2 = c(maxOct*cos(9*pi/24),maxOct*cos(10*pi/24), maxOct*cos(11*pi/24), maxOct*cos(13*pi/24), maxOct*cos(14*pi/24), maxOct*cos(15*pi/24)),
                            y1 = c(0,0,0,0,0,0), 
                            y2 = c(maxOct*sin(9*pi/24),maxOct*sin(10*pi/24), maxOct*sin(11*pi/24), maxOct*sin(13*pi/24), maxOct*sin(14*pi/24), maxOct*sin(15*pi/24)),
                            label = c(0.75,0.5,0.25,-0.25, -0.5, -0.75)
     )
     
     # maxOct*cos(pi/2)
     # maxOct*sin(pi/2)
     
     #construct dataframe to plot
     # xy <- data.frame(cbind(oct, lodes))
     # xyMean <- data.frame(cbind(meanOct, meanLode))
     
     #The Plot
     HsuPlot <- ggplot() + 
          coord_fixed(ratio = 1, ylim = c(-0.01,maxOct + 0.02), xlim = c(-maxOct/2 -0.01,maxOct/2 + 0.01)) +
          geom_arc(aes(x0=0,y0=0, r= maxOct, start = - pi/6, end = pi/6), linetype = 1, show.legend = FALSE) + 
          geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), linetype = 1, data = lines, show.legend = FALSE) +
          geom_text(aes(x=x2, y=y2,label = label, vjust = -0.5, family = "Helvetica"), size = 2.75, data = lines) +
          theme_void() +
          theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
     #plot the arcs and label them
     
     
     
     if (plotArcs == TRUE) {
          HsuPlot = HsuPlot + geom_arc(aes(x0 = 0, y0 = 0, r = r, start = start, end = end, alpha = 1), linetype = 1, color = "dark grey", data = arcs, show.legend = FALSE) + 
               geom_text(aes(x=textX, y=textY,label = r,angle = -30), size = 2.75, family = "Helvetica", nudge_x = 0.01, nudge_y = -0.01, data = arcs)
     }
     
     
     if (plotLines == TRUE) {
          #plot lines of equal lodes parameter and label them
          HsuPlot = HsuPlot + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color = "dark grey",linetype = 1, data = subLines, show.legend = FALSE) +
               geom_text(aes(x=x2, y=y2,label = label, vjust = -0.5, family = "Helvetica"),  size = 2.75, data = subLines) 
     }
     
     
     #geom_text(aes(x=x, y=y, label = labels, angle = angle, family = "Times"), size = 2.75, data = axisLabels)+
     #plot data
     # geom_point(aes(x=lodes, y=oct), data = xy, size = 1)+
     # geom_point(aes(x=meanLode, y=meanOct), data = xyMean, size = 4, color = "blue")+
     
     
     return(HsuPlot)
} 




##### This is a fully functioning Hsu-Nadai plot. It's main data input is a list of log-axis lengths of a set of ellipsoids. It computes the octahedral shear strain and lodes parameters and plots the points appropriately. 

#' @param plot a baseplot upon which to plot the dataset. E.g., could be a Hsu-Nadai plot with other data on it. If left as NULL, the data will be plotted on an auto-scaled blank HsuNadai plot.
#' @param logAs a list of vectors, each of length 3. The log-axis lengths of ellipsoids. 
#' @param maxOct the maximum octahedral shear strain in the plot. Only used if plot = NULL. If maxOct = NULL, the plot will automatically pick the nearest integer above the max octahedral shear strain in the dataset. 
#' @param alpha transparency. A vector of numbers between 0 and 1, either length 1 or the length of the dataset. If the alpha vector is not the length of the dataset, the first value will be repeated for the length of the dataset (e.g., if you input alpha = 1, then all points on the Hsu plot will have an alpha of 1). 
#' @param color the color of the points. A vector of strings, either R colors or hexidecmal, of length 1 or the length of the dataset. If the color vector is not the length of the dataset, the first value will be repeated for the length of the dataset (e.g., if you input color = "black", then all points on the Hsu plot will be black).
#' @param shape the shape of the points. See ggplot2 documentation for options. A vector of numbers of length 1 or the length of the dataset. If the shape vector is not the length of the dataset, the first value will be repeated for the length of the dataset (e.g., if you input shape = c(20), then all points on the Hsu plot will be filled circles).
#' @param size size of the points. vector of length 1 or length of dataset. follows same conventions as above.
#' @param annotateWith list of strings. Vector the length of the dataset. If NULL, no text annotations will be drawn. examples of annotations might be station names.  
#' @param hidePoints boolean. If TRUE, points will not be drawn. If annotateWith is not NULL, then text annotations will be centered on the not-drawn points and will have white frames around them, so that the annotations basically are the datapoints. 
plotHsuNadai = function(plot = NULL, logAsList = list(), maxOct = NULL, alpha = c(1), color = c("black"), shape = c(20), size = 1, annotateWith = NULL, hidePoints = FALSE, colorBy = NULL) {
     
     if (length(logAsList) > 1) {
          if (length(shape) != length(logAsList)) {
               shape = rep(shape[1],length(logAsList))
          }
          
          if (length(alpha) != length(logAsList)) {
               alpha = rep(alpha[1],length(logAsList))
          }
          
          if (length(color) != length(logAsList)) {
               color = rep(color[1],length(logAsList))
          }
          if (length(size) != length(logAsList)) {
               size = rep(size[1],length(logAsList))
          }
     }
     
     x <- function(logA) {-sqrt(1.5) * max(logA - sum(logA) / 3) - sqrt(1.5) * min(logA - sum(logA) / 3)}
     y <- function(logA) {sqrt(0.5) * max(logA - sum(logA) / 3) - sqrt(0.5) * min(logA - sum(logA) / 3)}
     

     
     if (is.null(plot) == TRUE){
          if (is.null(maxOct) == TRUE){
               autoMaxOct = 0
               for(i in 1:length(logAsList)){
                    
                    octs = sapply(logAsList[[i]],y)
                    if(autoMaxOct >= max(octs)){
                         autoMaxOct = autoMaxOct
                    } else {
                         autoMaxOct = max(octs)
                    }
               }
               plot = plotBlankHsuNadaiPlot(maxOct = ceiling(autoMaxOct))
          } else {
               plot = plotBlankHsuNadaiPlot(maxOct = maxOct)
          }
          
     }
     
     
     
     
     for(i in 1:length(logAsList)){
          logAs = logAsList[[i]]
          
          octs = sapply(logAs,y)
          lodes = sapply(logAs,x)
          xys =  data.frame(cbind(octs, lodes))
          

          
          
          
          if(hidePoints == FALSE) {
               if(is.null(colorBy) == TRUE){
                    for (n in 1:length(logAs)) {
                         logA = logAs[[n]]
                         
                         lode <- x(logA)
                         oct <- y(logA)
                         xy <- data.frame(cbind(oct, lode))
                         
                         plot = plot + geom_point(aes(x=lode, y=oct), data = xy, size = size[i], alpha = alpha[i], color = color[i], shape = shape[i])
                         
                         
                    }
               } else {
                    plot = plot + geom_point(aes(x=lodes, y=octs, color = colorBy), data = xys, size = size[1], alpha = alpha[1], shape = shape[1])
               }
               
               
               if (is.null(annotateWith) == FALSE) {
                    plot = plot + geom_text(data = xys, aes(x = lodes, y = octs, label = annotateWith), nudge_x = 0.1, nudge_y = 0.01, size = 2)
               }
          } else {
               if (is.null(annotateWith) == FALSE) {
                    plot = plot + geom_label(data = xys, aes(x = lodes, y = octs, label = annotateWith), size = 2)
               }
          }
          
     }
     
     
     
     
     return (plot)
     
}

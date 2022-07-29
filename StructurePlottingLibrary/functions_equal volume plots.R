#takes triangles in list format that the rgl plot3d uses and converts them into the matrix format that plotly mesh3d uses
trianglesRGLtoPLOTLY = function(trianglesRGL) {
     x = c()
     y = c()
     z = c()
     i = c()
     j = c()
     k = c()
     count = 0
     for (n in 1: length(trianglesRGL)) {
          t = trianglesRGL[[n]]
          x = c(x, t[[1]][1], t[[2]][1], t[[3]][1])
          y = c(y, t[[1]][2], t[[2]][2], t[[3]][2])
          z = c(z, t[[1]][3], t[[2]][3], t[[3]][3])
          i = c(i, count)
          j = c(j, count + 1)
          k = c(k, count + 2)
          
          count = count + 3
     }
     trainglesPLOTLY = data.frame(x,y,z,i,j,k)
}

plotEqualVolumeBlank = function(sphereOpacity = 0.05, sphereColor = "#efefef", axisColor = c("black", "black","black")) {
     # draw a sphere around the plot
     sphere = rayTetrahedralSphere(5)
     # define triangles from the sphere
     trianglesRGL = lapply(sphere, function(tri) lapply(tri, function(S) S * rotEqualVolumeRadius))
     # convert list of triangles to an x,y,z,i,j,k data frame
     triangles = trianglesRGLtoPLOTLY(trianglesRGL)
     
     # pre-define the plot axes to be blank
     axes <- list(title = "",
                  showgrid = FALSE,
                  zeroline = FALSE,
                  showline = FALSE,
                  showticklabels = FALSE)
     
     # pre-define east north up axes whose origin is at the center of the ball (0,0,0)
     eastAxis = data.frame("x" = c(rotEqualVolumeRadius,-rotEqualVolumeRadius), "y" = c(0,0), "z" = c(0,0))
     northAxis = data.frame("x" = c(0,0), "y" = c(rotEqualVolumeRadius,-rotEqualVolumeRadius), "z" = c(0,0))
     upAxis = data.frame("x" = c(0,0), "y" = c(0,0), "z" = c(rotEqualVolumeRadius,-rotEqualVolumeRadius))
     
     # pre-define the 3D annotations for east, north, and up
     eastAnn <- list(x = rotEqualVolumeRadius + 0.15,
                     y = 0,
                     z = 0,
                     xanchor = 'middle',
                     text = "East",
                     font = list(family = 'Helvetica',size = 16, color = 'rgba(0,0,0,1)'),
                     showarrow = FALSE)
     
     
     northAnn <- list(x = 0, 
                      y = rotEqualVolumeRadius + 0.15, 
                      z = 0, 
                      text = "North", 
                      xanchor = "middle", 
                      font = list(family = 'Helvetica',size = 16, color = 'rgba(0,0,0,1)'),
                      showarrow = FALSE)
     
     upAnn <- list( x = 0,
                    y = 0, 
                    z = rotEqualVolumeRadius + 0.1,
                    zanchor = 'bottom',
                    text = "Up",
                    font = list(family = 'Helvetica',size = 16,color = 'rgba(0,0,0,1)'),
                    showarrow = FALSE)
     
     
     #The plot
p <- plot_ly() %>%
          add_trace(data = triangles, x = ~x, y = ~y, z = ~z, i = ~i, j = ~j, k = ~k,
             type = "mesh3d", 
             facecolor = rep(I(sphereColor), nrow(triangles)), 
             opacity = sphereOpacity, 
             lighting = list(specular = 2), 
             lightposition = list(x = 10, y = 10, z = 10))  %>%
          # add_axes
          
          add_trace(data = eastAxis, name = "East", x= ~x, y = ~y, z = ~z, type = "scatter3d", mode = "lines+markers", marker = list(color = I(axisColor[1]), size = 4), line = list(color = I(axisColor[1]),width = 5), opacity = 1,
                    inherit = FALSE) %>%
          add_trace(data = northAxis, name = "North", x= ~x, y = ~y, z = ~z, type = "scatter3d",mode = "lines+markers", marker = list(color = I(axisColor[2]), size = 4), line = list(color = I(axisColor[2]),width = 5), opacity = 1,
                    inherit = FALSE) %>%
          add_trace(data = upAxis, name = "Up", x= ~x, y = ~y, z = ~z, type = "scatter3d",mode = "lines+markers", marker = list(color = I(axisColor[3]), size = 4), line = list(color = I(axisColor[3]),width = 5), opacity = 1,
                    inherit = FALSE) %>% layout(showlegend = FALSE) %>%
          layout(
               scene = list(
                    xaxis = axes,
                    yaxis = axes,
                    zaxis = axes,
                    aspectratio = list(
                         x = 1,
                         y = 1,
                         z = 1
                    ),
                    # camera = list(
                    #      center = list(x = 0,
                    #           y = 0,
                    #           z = 0
                    #      ),
                    #      eye = list(
                    #           x = 1,
                    #           y = -1,
                    #           z = 1
                    #      ),
                    #      up = list(
                    #           x = 0,
                    #           y = 0,
                    #           z = 1
                    #      )
                    # ),
                    annotations = list(
                         eastAnn, northAnn, upAnn
                    ))
          )
          p
          
}



plotEqualVolumeRotations = function(points = NULL, trianglesRot = NULL, group = oriLineInPlaneGroup, triColor = list(I("#999999")), color = list(I("black")), sphereOpacity = 0.05, sphereColor = "#efefef", axisColor = c("black", "black","black"), pointSize = 2) {
     # account for symmetry of group
     if (is.null(points) == FALSE) {
         equalVolumePoints = list()
        for (i in 1:length(points)){
            points[[i]] <- oriSymmetrizedRotations(points[[i]], group)
            equalVolumePoints[[i]] = data.frame(t(sapply(points[[i]], function(s) rotEqualVolumeFromMatrix(s))))
        }
     }
    for (i in 1:length(color)){
        if (length(color[[i]]) > 1) {
            color[[i]] = rep(color[[i]], length(group))
        }
    }    

    for (i in 1:length(triColor)){
        if (length(triColor[[i]]) > 1) {
            triColor[[i]] = rep(triColor[[i]], length(group))
        }
    }    
     
     if (is.null(trianglesRot) == FALSE) {
         triangles = list()
         f <- function(g, ps) {lapply(ps, function(p) g %*% p)}
         for (i in 1:length(trianglesRot)){
             trianglesSymm <- Reduce(c, lapply(group, function(g) lapply(trianglesRot[[i]], function(ps) f(g, ps))))
             trianglesEV <- lapply(trianglesSymm, function(tri) {lapply(tri, rotEqualVolumeFromMatrix)})
             triangles[[i]] = trianglesRGLtoPLOTLY(trianglesEV) 
         }
     }
    
     
     
     p = plotEqualVolumeBlank(sphereOpacity = sphereOpacity, sphereColor = sphereColor, axisColor = axisColor)
     
     if (is.null(points) == FALSE & is.null(trianglesRot) == FALSE) {
     
         for (i in 1:length(equalVolumePoints)) {
             p = p %>%
             add_markers(data = equalVolumePoints[[i]], name = paste("points",i) ,
                         x= ~X1, y = ~X2, z = ~X3, type = "scatter3d",
                         #opacity = 1, 
                         color = color[[i]], 
                         marker = list(size = pointSize),
                         inherit = FALSE)
         } 
         for (i in 1:length(triangles)){
         p = p %>%
           add_trace(data = triangles[[i]], x = ~x, y = ~y, z = ~z, i = ~i, j = ~j, k = ~k,
               type = "mesh3d", 
               facecolor = rep(triColor[[i]], nrow(triangles[[i]])), 
               opacity = 1, 
               lighting = list(specular = .5), 
               lightposition = list(x = 10, y = 10, z = 10))
         }
     }
     
     if (is.null(points) == FALSE & is.null(trianglesRot) == TRUE) {
         
         for (i in 1:length(equalVolumePoints)) {
             p = p %>%
             add_markers(data = equalVolumePoints[[i]], name = paste("points",i) ,
                         x= ~X1, y = ~X2, z = ~X3, type = "scatter3d",
                         #opacity = 1, 
                         color = color[[i]], 
                         marker = list(size = pointSize),
                         inherit = FALSE)
         } 
     }
     
     p
          
}

plotEqualVolumeEllipsoidRotations = function(ellipsoidRotations, trianglesRot = NULL, group = oriLineInPlaneGroup, triColor = list(I("#999999")), color = list(I("black")), sphereOpacity = 0.05, sphereColor = "#efefef", axisColor = c("black", "black","black"), pointSize = 2) {
    
    # reorder rows-- row1 = pole to foliation; row2 = lineation, row3 = third orthogonal
    rotations = lapply(ellipsoidRotations, function(rot){
        row1 = rot[3,]
        row2 = rot[1,]
        row3 = rot[2,]
        rotation = rbind(row1,row2,row3)
    })
 
    
    plotEqualVolumeRotations(list(rotations),  trianglesRot = trianglesRot, group = group, triColor = triColor, color = color, sphereOpacity = sphereOpacity, sphereColor = sphereColor, axisColor = axisColor, pointSize = pointSize)   
}








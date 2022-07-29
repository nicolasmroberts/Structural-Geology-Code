# This tutorial shows you how to use the plotting functions I wrote built on the ggplot2 package. The functions are built produce equal area nets that in a highly flexible, customizable progamatic way. Using the underlying ggplot2 engine, you can interatively built up your plots with different layers, color points by independent parameters or customize the colors by hand.

# Import some dependent libraries. If you don't have all the packages needed, you'll need to install them before running these lines of code. Open up the two files listed below, and R should tell you whether you need to install anything.
source("GeologyGeometryLibrary/all.R")
source("StructurePlottingLibrary/all.R")


# Import example dataset
follins = read.csv("exampleFollins.csv")


##### Quick plot
# This dataset has easting, northing, strike, dip, rake, and feature type (foliation, vein, bedding). First, let's just worry about the strikes and dips. 
# You can directly plot strikes and dips using the function plotEqualAreaNetSTDP. 

# As an input, it needs a list of strike-dip pairs and/or trend-plunge pairs. This line of code produces such a list of strikes and dips.
SDs = lapply(1:length(follins$strike), function(row) c(follins$strike[row], follins$dip[row]))

# Now plot. 
plotEqualAreaNetSDTP(sdPairs = SDs)





##### More indepth
##### 
# The above function is really a shortcut to the main function, and has pretty much no options besides changing the color. To access the full function, you need to change your data from degrees (strikes and dips) to 3D cartesian vectors. The following line of code appends a column onto the follins dataframe of cartesian vectors that represent the pole to foliation. 
follins$poles = lapply(1:length(follins$strike), function(row) geoCartesianFromStrikeDipDeg(c(follins$strike[row], follins$dip[row])))

# Before we get going, let's look at some options for a blank equal area net.
# The Default blank net.
plotBlankEqualAreaNet()

# Here are all the options you have. Play around with these parameters. 
plotBlankEqualAreaNet(centerCross = TRUE, cardinalTicks = TRUE, cardinalLabels = TRUE, circleLineWeight = 3, fillColor = "light blue", lineColor = "green")

#


# Now, let's plot those poles to foliation. 
plotEqualAreaNet(points = list(follins$poles))


# We can do longs of things with this. Let's change the points from squares to circles. The shape number corresponds to GGPLOT2's shape numbers (google them!)
plotEqualAreaNet(points = list(follins$poles),
                 shape = 21)

# If you want to change the look of the equal area net, then do the following:
plotEqualAreaNet(plot = plotBlankEqualAreaNet(centerCross = FALSE), # This line can contain ANY GGPLOT object.
                 points = list(follins$poles),
                 shape = 21)

# But maybe we want to plot poles to foliation and lineation on the sample plot! let's load up the lineations as cartesian vectors
# compute a matrix (first row: lineation, second row: pole to plane, third row: strike vector)
follins$rotation = lapply(1:length(follins$strike), function(row) geoCartesianFromStrikeDipRakeDeg(c(follins$strike[row], follins$dip[row], follins$rake[row])))

# extract just the first row (the lineation)
follins$lines = lapply(follins$rotation, function(rot) rot[1,])

# And now plot both the poles to foliation and the lineations. To do that, just list both under the points input.
plotEqualAreaNet(points = list(follins$poles,
                               follins$lines))

# You can choose different shapes for each point dataset. 
plotEqualAreaNet(points = list(follins$poles,
                               follins$lines),
                 shape = c(21,22),
                 color = c("blue","orange"),
                 edgeColor = c("white", "black"),
                 pointAlpha = c(0.5,1))

# You might want to plot the poles as great circles instead of points!
plotEqualAreaNet(points = list(follins$lines),
                 curves = list(follins$poles))


# You might want to plot the poles as great circles instead of points!
plotEqualAreaNet(points = list(follins$lines),
                 curves = list(follins$poles),
                 color = "orange",
                 edgeColor = "black",
                 pointAlpha = 1,
                 curveColor = "blue",
                 curveAlpha = 0.5)

#Another think you might like to do is color the points by some parameter. Say, their struture time (bedding or foliation)
plotEqualAreaNet(points = list(follins$lines),
                 curves = list(follins$poles),
                 colorBy = follins$featureType,
                 showLegend = TRUE)

# Or color by easting
plotEqualAreaNet(points = list(follins$lines),
                 curves = list(follins$poles),
                 colorBy = follins$easting,
                 showLegend = TRUE)


# Or color by northing (this time, poles and lines as points)
plotEqualAreaNet(points = list(follins$lines,
                               follins$poles),
                 shape = c(22,21),
                 colorBy = follins$easting,
                 showLegend = TRUE)


##### Kamb Contouring
# To generate contours, use this function. Breaks is the sigma level you want to draw (change as you see fit). numNonAdapt is how detailed (or pixelated) the contour regions will be. Change as you see fit (but higher numbers will be much slower)
plotKambContours(points = follins$poles, numNonAdapt = 4, breaks = c(1,3,6)) 


# Add points. You do this by using the kamb contour plot from above as the baseplot in the function.
plotEqualAreaNet( plot = plotKambContours(points = follins$poles, numNonAdapt = 5, breaks = c(0,3,9,22)) ,
                  points = list(follins$poles,
                                follins$lines),
                  shape = c(21,22),
                  color = c("blue","orange"),
                  edgeColor = c("white", "black"),
                  pointAlpha = c(0.5,1))

# Change the color scale for contours. 
plotEqualAreaNet( plot = plotKambContours(points = follins$poles, numNonAdapt = 5, breaks = c(0,3,9,22)) ,
                  points = list(follins$poles,
                                follins$lines),
                  shape = c(21,22),
                  color = c("blue","orange"),
                  edgeColor = c("white", "black"),
                  pointAlpha = c(0.5,1)) + scale_color_viridis_d() + scale_fill_viridis_d()

# You can accomplish the same thing by building the plot in a few lines. 
kambPlot = plotKambContours(points = follins$poles, numNonAdapt = 5, breaks = c(0,3,9,22))

plotEqualAreaNet( plot = kambPlot,
                  points = list(follins$poles,
                                follins$lines),
                  shape = c(21,22),
                  color = c("blue","orange"),
                  edgeColor = c("white", "black"),
                  pointAlpha = c(0.5,1)) + scale_color_viridis_d() + scale_fill_viridis_d()




##### Example with bootstrapping

# In this example, we will bootstrap the lineation data, get the mean and confidence region, and plot all of it up in one spot.
lineBoots = lineBootstrapInference(follins$lines,numBoots = 1000, numPoints = 100)

# Plot the bootstraps
plotEqualAreaNet(points = list(lineBoots$us))

# That doesn't looks so good, because all the points are on top of each other, and the white border of the points makes the edges disappear. Let's fix that. Also let's make the points be almost complete transparent so we can get an idea of the density of the bootstrap point cloud. 
plotEqualAreaNet(points = list(lineBoots$us),
                 shape = c(21), edgeColor = "black",
                 pointAlpha = 0.002)

# Draw the 95% confidence region for the population mean
plotEqualAreaCI(polygonPoints = list(lineBoots$points), edgeColor = "black", fillColor = NA)

#And put it all together
layer1 = plotEqualAreaNet(points = list(lineBoots$us),
                           shape = c(21), edgeColor = "black",
                           pointAlpha = 0.002)

layer2 = plotEqualAreaNet(plot = layer1, points = list(follins$lines),
                          pointAlpha = 0.2,
                          edgeColor = "black")

layer3 = plotEqualAreaCI(plot = layer2,polygonPoints = list(lineBoots$points), edgeColor = "white", fillColor = NA )

layer3




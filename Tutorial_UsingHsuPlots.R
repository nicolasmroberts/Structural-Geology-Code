# This tutorial shows you how to use the plotting functions I wrote built on the ggplot2 package. The functions are built produce Hsu-Nadai plots that are highly flexible, customizable and progamatic. Using the underlying ggplot2 engine, you can interatively built up your plots with different layers, color points by independent parameters or customize the colors by hand.

# Import some dependent libraries. If you don't have all the packages needed, you'll need to install them before running these lines of code. Open up the two files listed below, and R should tell you whether you need to install anything.
source("GeologyGeometryLibrary/all.R")
source("StructurePlottingLibrary/all.R")

##### Import data from Roberts and Tikoff, 2020, Precambrian research
#load csv file
dataFrame <- geoDataFromFile("example_amsData.csv", separator = ",")

#This computes useful measures from the raw AGICO csv. This is the main variable which all the statistics will draw from. 
amsData <- geoEllipsoidDataFromAGICOFile("example_amsData.csv", separator = ",", sapply(dataFrame$Km, function(s) s*10^(-6)), doNormalize=TRUE)

#Loads a csv file with the location (easting, northing) information. The table needs to be formatted with headings "name", "easting", and "northing"--all lower case, no spaces. The table can include more stations than the amsData have, but must inlcude all stations in amsData. More details below. 
locData <- read.table("example_sampleLocations.csv", sep = ",", header=TRUE,stringsAsFactors = FALSE)

nameLength <- 9

#this calls the locationPairing function and apppends location data to the amsDataset. 
amsData <- locationPairing(amsData, locData, nameLength)
amsData$easting = as.numeric(amsData$easting) # make sure the coordinates are numerics (they default to character class)
amsData$northing = as.numeric(amsData$northing)  # make sure the coordinates are numerics (they default to character class)

# Remove duplicate station names
stationList <- unique(amsData$station, incomparables = FALSE)
amsDataMeans = computeStationMeans(amsData, stationList, isSimpleFeature = FALSE)



##### Hsu-Nadai Plot
# You can plot a blank HsuNadai plot
plotBlankHsuNadaiPlot()

# Here are all the options for the plot
plotBlankHsuNadaiPlot(plotArcs = TRUE, numArcs = 1, plotLines = FALSE, maxOct = 2)

# Here is a basic plot. The only required input is a list of LogAs (log-axis lengths of the ellipsoids)
plotHsuNadai(logAs = list(amsDataMeans$logA))

# The plot automatically scales the radial axis to the nearest integer above the largest point. If you want to tweak that, you can use maxOct.
plotHsuNadai(logAs = list(amsDataMeans$logA), maxOct = 0.25)

# You can change the size of the points
plotHsuNadai(logAs = list(amsDataMeans$logA), maxOct = 0.25, size = 3)

# You can also layer multiple datasets.
plotHsuNadai(logAsList = list(amsDataMeans$logA,
                              amsData$logA), 
             maxOct = 0.25,
             size = c(3,1), 
             alpha = c(1,0.3))

# Or you can color the points by some parameter.
plotHsuNadai(logAs = list(amsDataMeans$logA), 
             maxOct = 0.25, 
             size = 3,
             colorBy = amsDataMeans$Km)

# And because it is a ggplot object, you can always change the colorscale
plotHsuNadai(logAs = list(amsDataMeans$logA), 
             maxOct = 0.25, 
             size = 3,
             colorBy = log10(amsDataMeans$Km)) + scale_color_viridis()

# Or color the points by northing
plotHsuNadai(logAs = list(amsDataMeans$logA), 
             maxOct = 0.25, 
             size = 3,
             colorBy = amsDataMeans$northing) + scale_color_viridis()

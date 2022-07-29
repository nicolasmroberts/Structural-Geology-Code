


downDipFollin = function(rotation) {
     strikeDipRakeDeg = geoStrikeDipRakeDegFromRotation(rotation)
     downDipTrendPlungeDeg = c(mod(strikeDipRakeDeg[1] + 90, 360), strikeDipRakeDeg[2])
     downDipTrendPlungeRotation = geoRotationFromStrikeDipTrendPlungeDeg(c(strikeDipRakeDeg[1], strikeDipRakeDeg[2], downDipTrendPlungeDeg))
     downDipTrendPlungeRotation
}

distanceFromDownDip = function(rotation) {
     downDipFollin = downDipFollin(rotation)
     line = downDipFollin[2,]
     distance = lineDistance(line, rotation[2,]) / degree
     strikeDipRakeDeg = geoStrikeDipRakeDegFromRotation(rotation)
     if (strikeDipRakeDeg[3] > 90) {
          distance = -distance
     }
     distance
}

rotationFromTwoRotations = function(rotation1, rotation2) {
        rotation = t(rotation2) %*% rotation1
}

bootStrappedRotationDifference = function(bootStraps1, bootStraps2) {
        rotations = lapply(1:length(bootStraps1), function(s) rotationFromTwoRotations(bootStraps1[[s]], bootStraps2[[s]]))
        rotations
}


limitRot95 = function(rotBoots) {
mahalanobisDs = lapply(1: length(rotBoots$bootstraps), function(s) rotMahalanobisNorm(rotBoots$bootstraps[[s]], rotBoots$center, rotBoots$leftCovarInv) )
rot95 = list()
counter = 1
for (i in 1:length(rotBoots$bootstraps)){
        if (mahalanobisDs[[i]] < rotBoots$q095) {
                rot95[[counter]] = rotBoots$bootstraps[[i]]
                counter = counter + 1
        }
}
rot95
}



#data frame from readCSV and a critirion for which rows contain foliation-lineation pairs
rotationsFromeCSVframe = function(dataFrame, FollinCrit) {
     
     Follins <- data.frame(easting = dataFrame$easting[FollinCrit], 
                           northing = dataFrame$northing[FollinCrit], 
                           strike = dataFrame$strike[FollinCrit],
                           dip = dataFrame$dip[FollinCrit], 
                           rake = dataFrame$rake[FollinCrit])
     
     rotations = list()
     
     for (i in 1:length(Follins$easting)){
          strike = Follins$strike[[i]]
          dip = Follins$dip[[i]]
          rake = Follins$rake[[i]]
          
          rotations[[i]] <- geoRotationFromStrikeDipRakeDeg(c(strike, dip, rake))
     }
     
     return(rotations)
}






rotationsFromeCSVframeFolds = function(dataFrame, Crit) {
     
     Folds <- data.frame(easting = dataFrame$easting[Crit], 
                         northing = dataFrame$northing[Crit], 
                         strike = dataFrame$axialPlaneStrike[Crit],
                         dip = dataFrame$axialPlaneDip[Crit], 
                         trend = dataFrame$hingeTrend[Crit],
                         plunge = dataFrame$hingePlunge[Crit])
     
     rotations = list()
     
     for (i in 1:length(Folds$easting)){
          strike = Folds$strike[[i]]
          dip = Folds$dip[[i]]
          trend = Folds$rake[[i]]
          plunge = Folds$plunge[[i]]
          
          rotations[[i]] <- geoRotationFromStrikeDipTrendPlungeDeg(c(strike, dip, trend,plunge))
     }
     
     return(rotations)
}



cartesianFromeCSVframeStrikeDip = function(dataFrame, Crit) {
     
     Follins <- data.frame(easting = dataFrame$easting[Crit], 
                           northing = dataFrame$northing[Crit], 
                           strike = dataFrame$strike[Crit],
                           dip = dataFrame$dip[Crit] )
     
     cartesianVecs = list()
     
     for (i in 1:length(Follins$easting)){
          strike = Follins$strike[[i]]
          dip = Follins$dip[[i]]
          
          cartesianVecs[[i]] <- geoCartesianFromStrikeDipDeg(c(strike, dip))
     }
     
     return(cartesianVecs)
}




#one Regression does a regression based on an azimuth (e.g. how well does northeasting (045) explain variation in my data?) This kernal geodetic regression is analogous to a best-fit line, and thus describes a steady variation, not a sudden change in orientation.
oneRegression <-
        function(follins,
                 domain,
                 pValuePerms,
                 directionality) {
                regression <-
                        oriGeodesicRegression(
                                cos((directionality) * degree) * follins$northing[domain] + sin((directionality) *
                                                                                                          degree) * follins$easting[domain],
                                follins$rotation[domain],
                                oriLineInPlaneGroup,
                                numSteps = 10000
                        )
                
                if (pValuePerms > 0) {
                        RSquareds <-
                                oriGeodesicRegressionPermutations(
                                        cos((directionality) * degree) * follins$northing[domain] + sin((directionality) *
                                                                                                                  degree) * follins$easting[domain],
                                        follins$rotation[domain],
                                        numPerms = pValuePerms,
                                        group = oriLineInPlaneGroup
                                )
                        sum(RSquareds > regression$rSquared)
                        p <- sum(RSquareds > regression$rSquared) / length(RSquareds)
                } else {
                        p <- 0
                }
                regressionStats <- zeros(1, 5)
                regressionStats[1, 1] = directionality
                regressionStats[1, 2] = regression$error
                regressionStats[1, 3] = regression$minEigenvalue
                regressionStats[1, 4] = regression$rSquared
                regressionStats[1, 5] = p
                names(regressionStats) = c("azimuth", "error", "minEigenValue", "R^2", "P")
                regression$stats <- regressionStats
                return(regression)
        }

#regressionSweep does a series of directional regressions from 0-180°, given an increment, and also calculates the pValue for each regression. BEWARE: doing this with p-values can take days-weeks-or-months of computing time. I suggest first running the function with degreeIncrement = 10 and pValuePerms=0. This may still take a couple hours, but will give you a first order picture of whether it is interesting to proceed.
regressionSweep <- function(follins, degreeIncrement = 10, pValuePerms = 0) {
                intervals <- 180 / degreeIncrement
                regressionDeg = seq(0,180, degreeIncrement)
                # Geodesic regression of pole vs. northing in domain 4. Check that error is 0 and minEigenvalue is positive.
                
                reg= oriRescaledGeodesicRegression(cos(regressionDeg[1] * degree) * follins$northing + sin(regressionDeg[1] * degree) * follins$easting,
                                                   follins$rotation,  oriLineInPlaneGroup, numSteps = 10000)
                
                regressionData = data.frame("azimuth" = regressionDeg[1], "rsquared" = reg$rSquared, "error" = reg$error, "minEig" = reg$minEigenvalue)
                
                
                for ( i in 2: length(regressionDeg)) {
                        reg = oriRescaledGeodesicRegression(cos(regressionDeg[i] * degree) * follins$northing + sin(regressionDeg[i] * degree) * follins$easting,
                                                           follins$rotation,  oriLineInPlaneGroup, numSteps = 10000)
                        
                        regFrame = data.frame("azimuth" = regressionDeg[i], "rsquared" = reg$rSquared, "error" = reg$error, "minEig" = reg$minEigenvalue)
                        regressionData = rbind(regressionData, regFrame)
                }
        
                regressionData
                        
                #         if (pValuePerms > 0) {
                #                 RSquareds <-
                #                         oriGeodesicRegressionPermutations(
                #                                 cos((v[i, 1] - 1) * degreeIncrement * degree) * follins$northing[domain] + sin((v[i, 1] -1) * degreeIncrement * degree) * follins$easting[domain],
                #                                 follins$rotation[domain],
                #                                 numPerms = pValuePerms,
                #                                 group = oriLineInPlaneGroup)
                #                 length(RSquareds)
                #                 sum(RSquareds > regressionTemp$rSquared)
                #                 p <-
                #                         sum(RSquareds > regressionTemp$rSquared) / length(RSquareds)
                #         }
                #         else {
                #                 p <- "nan"
                #         }
                #         
                #         regressionStats[i, 5] = p
                # }
                # regressions[[i + 1]] <- regressionStats
                # return(regressions)
        }










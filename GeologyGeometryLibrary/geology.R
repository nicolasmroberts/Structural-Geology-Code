


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This is the part of the library that is designed specifically for geologic applications. In short, our strategy is to convert geologic data into abstract math, do computations there, and then convert the answers back into a form that a geologist can understand. For example, a geologist might think of a fault plane in terms of strike and dip (in degrees), but we compute with its pole vector expressed in Cartesian coordinates. By default, our x-y-z coordinates point east, north, up. Any function or variable that deals with degrees has 'Deg' in its name, like a skull and crossbones on a bottle of poison.



### COORDINATE TRANSFORMATIONS ###

geoSphericalFromTrendPlungeDeg <- function(tpDeg) {
	rho <- 1.0
	phi <- tpDeg[2] * degree + pi / 2.0
	theta <- (pi / 2.0 - tpDeg[1] * degree) %% (2.0 * pi)
	c(rho, phi, theta)
}

geoTrendPlungeDegFromSpherical <- function(rpt) {
	trendDeg <- (90.0 - rpt[3] / degree) %% 360.0
	plungeDeg <- rpt[2] / degree - 90.0
	c(trendDeg, plungeDeg)
}

#' Cartesian coordinates from trend and plunge in degrees.
#' 
#' Inverse to geoTrendPlungeDegFromCartesian, up to periodicity of angles.
#' @param tpDeg A 2D real vector. The trend and plunge of a ray or line, in degrees.
#' @return A 3D real vector (unit).
geoCartesianFromTrendPlungeDeg <- function(tpDeg) {
	cartesianFromSpherical(geoSphericalFromTrendPlungeDeg(tpDeg))
}

#' Trend and plunge in degrees, from Cartesian coordinates.
#' 
#' Inverse to geoCartesianFromTrendPlungeDeg, up to periodicity of angles.
#' @param xyz A 3D real vector (unit). A ray or line.
#' @return A 2D real vector. Trend and plunge, in degrees.
geoTrendPlungeDegFromCartesian <- function(xyz) {
	geoTrendPlungeDegFromSpherical(sphericalFromCartesian(xyz))
}

geoSphericalFromStrikeDipDeg <- function(sdDeg) {
	rho <- 1.0
	phi <- sdDeg[2] * degree
	theta <- 2.0 * pi - sdDeg[1] * degree
	c(rho, phi, theta)
}

geoStrikeDipDegFromSpherical <- function(rpt) {
	if (rpt[2] <= pi / 2.0) {
		sDeg <- (360.0 - rpt[3] / degree) %% 360.0
		dDeg <- rpt[2] / degree
	} else {
		sDeg <- (360.0 - (rpt[3] + pi) / degree) %% 360.0
		dDeg <- (pi - rpt[2]) / degree
	}
	c(sDeg, dDeg)
}

#' Cartesian coordinates of a pole vector, from strike and dip in degrees.
#' 
#' Assumes the right-hand rule for strike and dip. Inverse to geoStrikeDipDegFromCartesian, up to periodicity of angles.
#' @param sdDeg A 2D real vector. The strike and dip of a plane, in degrees.
#' @return A 3D real vector (unit). The pole.
geoCartesianFromStrikeDipDeg <- function(sdDeg) {
	cartesianFromSpherical(geoSphericalFromStrikeDipDeg(sdDeg))
}

#' Strike and dip in degrees, from Cartesian coordinates of the pole vector.
#' 
#' Obeys the right-hand rule for strike and dip. Inverse to geoCartesianFromStrikeDipDeg, up to periodicity of angles.
#' @param xyz A 3D real vector (unit). The pole to the plane.
#' @return A 2D real vector. Strike and dip, in degrees.
geoStrikeDipDegFromCartesian <- function(xyz) {
	geoStrikeDipDegFromSpherical(sphericalFromCartesian(xyz))
}

# Cartesian coordinates for rake, pole, and strike, from strike, dip, and rake in degrees.
# Assumes the right-hand rule for strike and dip. Rake (also called pitch) is measured down-dip from the strike --- in other words, clockwise when viewed from above. Inverse to geoStrikeDipRakeDegFromCartesian, up to periodicity of angles. The return value is a real 3x3 matrix, usually not special orthogonal. The first row is the rake vector. The second row is the pole to the plane. The third row is the strike vector. The second row should be roughly perpendicular to the other two, but is often not exactly perpendicular due to imprecision in geologic data sets.
geoCartesianFromStrikeDipRakeDeg <- function(sdrDeg) {
  s <- -sdrDeg[1] * degree;
  s <- matrix(c(cos(s), sin(s), 0, -sin(s), cos(s), 0, 0, 0, 1), 3, 3)
  d <- sdrDeg[2] * degree;
  d <- matrix(c(cos(d), 0, -sin(d), 0, 1, 0, sin(d), 0, cos(d)), 3, 3)
  r <- -sdrDeg[3] * degree + pi / 2;
  r <- matrix(c(cos(r), sin(r), 0, -sin(r), cos(r), 0, 0, 0, 1), 3, 3)
  m <- t(s %*% d %*% r)
  rbind(m[1,], m[3,], s[,2])
}

# Strike, dip, and rake in degrees, from Cartesian coordinates for rake, pole, and strike.
# Inverse to geoCartesianFromStrikeDipRakeDeg, up to periodicity of angles.
geoStrikeDipRakeDegFromCartesian <- function(rakePoleStrike) {
  strikeDeg <- geoTrendPlungeDegFromCartesian(rakePoleStrike[3,])[1]
  dipDeg <- 90 - geoTrendPlungeDegFromCartesian(-rakePoleStrike[2,])[2]
  cosine <- cos(strikeDeg * degree)
  sine <- sin(strikeDeg * degree)
  s <- matrix(c(cosine, sine, 0, -sine, cosine, 0, 0, 0, 1), 3, 3)
  cosine <- cos(-dipDeg * degree)
  sine <- sin(-dipDeg * degree)
  d <- matrix(c(cosine, 0, -sine, 0, 1, 0, sine, 0, cosine), 3, 3)
  rakeDeg <- geoTrendPlungeDegFromCartesian(d %*% s %*% rakePoleStrike[1,])[1]
  c(strikeDeg, dipDeg, rakeDeg)
}

#' Orientation matrix from strike, dip, and rake in degrees.
#' 
#' Assumes the right-hand rule for strike and dip. Rake (also called pitch) is measured down-dip from the strike --- in other words, clockwise when viewed from above. Inverse to geoStrikeDipRakeDegFromRotation, up to periodicity of angles. When using this function to describe faults with slip, you often want to follow up with geoPoleVorticityFromPoleHanging.
#' @param sdrDeg A real 3D vector. The strike, dip, and rake, in degrees.
#' @return A real 3x3 matrix (special orthogonal). The first row is the pole to the plane. The second row is the ray in the plane described by the rake. (The third row is the cross product of the first two.)
geoRotationFromStrikeDipRakeDeg <- function(sdrDeg) {
  cart <- geoCartesianFromStrikeDipRakeDeg(sdrDeg)
  rotProjectedMatrix(rbind(cart[2,], cart[1,], cross(cart[2,], cart[1,])))
}

#' Strike, dip, and rake in degrees, from orientation matrix.
#' 
#' Assumes the right-hand rule for strike and dip. Rake (also called pitch) is measured down-dip from the strike --- in other words, clockwise when viewed from above. Inverse to geoRotationFromStrikeDipRakeDeg, up to periodicity of angles. When using this function to describe faults with slip, you often want to prepare with geoPoleHangingFromPoleVorticity.
#' @param rot A real 3x3 matrix (special orthogonal). The first row is the pole to the plane. The second row is a ray within the plane, whose rake is desired. (The third row is the cross product of the first two.)
#' @return A real 3D vector. The strike, dip, and rake, in degrees.
geoStrikeDipRakeDegFromRotation <- function(rot) {
  pole <- -lower(rot[1,])
  strike <- c(-pole[[2]], pole[[1]], 0)
  rake <- rot[2,]
  geoStrikeDipRakeDegFromCartesian(rbind(rake, pole, strike))
}

# Returns a tuple (strikeDeg, dipDeg, trendDeg, plungeDeg) describing the first two rows of the given rotation matrix. Inverse to geoRotationFromStrikeDipTrendPlungeDeg (ignoring angle conventions).
geoStrikeDipTrendPlungeDegFromRotation <- function(r) {
  planeDeg <- geoStrikeDipDegFromCartesian(r[1,])
  lineDeg <- geoTrendPlungeDegFromCartesian(r[2,])
  c(planeDeg, lineDeg)
}

# Returns a rotation matrix whose first row is the pole vector specified by the strike and dip in degrees, and whose second row is the vector specified by the trend and plunge in degrees.
geoRotationFromStrikeDipTrendPlungeDeg <- function(sdtpDeg) {
  plane <- geoCartesianFromStrikeDipDeg(sdtpDeg[1:2])
  line <- geoCartesianFromTrendPlungeDeg(sdtpDeg[3:4])
  other <- cross(plane, line)
  r <- t(matrix(c(plane, line, other), 3, 3))
  rotProjectedMatrix(r)
}

# Reports the axis and angle of the rotation. The axis is lower-hemisphere, so the angle may be negative.
geoTrendPlungeAngleDegFromMatrix <- function(r) {
  ua <- rotAxisAngleFromMatrix(r)
  if (ua[[3]] > 0)
    c(geoTrendPlungeDegFromCartesian(lower(ua[1:3])), -ua[[4]] / degree)
  else
    c(geoTrendPlungeDegFromCartesian(ua[1:3]), ua[[4]] / degree)
}

#' Sanitizing an orientation of a fault-with-slip.
#' 
#' The issue here is that geologists like describing fault slips in terms of the movement direction of the hanging wall, but the hanging wall is undefined for vertical faults. A much better-behaved convention is to describe the direction of slip using a vorticity vector. Inverse to geoPoleHangingFromPoleVorticity.
#' @param rot A real 3x3 matrix (special orthogonal). The first row is the pole to a fault plane. The second row is the movement direction of the hanging wall. (The third row is the cross product of the other two.)
#' @return A real 3x3 matrix (special orthogonal). The first row is the downward-pointing pole to the fault. The second row is the vorticity vector of the slip. (The third row is the cross product of the other two.)
geoPoleVorticityFromPoleHanging <- function(rot) {
  if (rot[[1, 3]] < 0)
    rbind(rot[1,], -rot[3,], rot[2,])
  else
    rbind(-rot[1,], rot[3,], rot[2,])
}

#' Desanitizing an orientation of a fault-with-slip.
#' 
#' Inverse to geoPoleVorticityFromPoleHanging.
#' @param rot A real 3x3 matrix (special orthogonal). The first row is the pole to a fault plane. The second row is the vorticity vector of the fault slip. (The third row is the cross product of the other two.)
#' @return A real 3x3 matrix (special orthogonal). The first row is the downward-pointing pole to the fault. The second row is the movement direction of the hanging wall. (The third row is the cross product of the other two.)
geoPoleHangingFromPoleVorticity <- function(rot) {
  if (rot[[1, 3]] < 0)
    rbind(rot[1,], rot[3,], -rot[2,])
  else
    rbind(-rot[1,], -rot[3,], -rot[2,])
}

# Positive (negative) longitudes are east (west). Positive (negative) latitudes are north (south).
# The positive x-axis is at the equator and prime meridian. The positive z-axis is at the north pole.
geoCartesianFromLatLongDeg <- function(latLongDeg) {
  cartesianFromSpherical(c(1, pi / 2 - latLongDeg[[1]] * degree, latLongDeg[[2]] * degree))
}

geoLatLongDegFromCartesian <- function(u) {
  sph <- sphericalFromCartesian(u)
  c(90 - sph[[2]] / degree, sph[[3]] / degree)
}

#' Average of points on the Earth's surface, in terms of latitude and longitude.
#' 
#' Assumes that the Earth is a perfect sphere. Positive latitudes are north, while negative latitudes are south. Positive longitudes are east, while negative longitudes are west. Handles the prime antimeridian (near international date line) and polar regions well.
#' @param latLongDegs A list of 2D real vectors, where each vector is a latitude-longitude pair in decimal degrees.
#' @return A 2D real vector, consisting of the latitude and longitude of the mean location, in decimal degrees.
geoMeanLatLongDeg <- function(latLongDegs) {
  carts <- lapply(latLongDegs, geoCartesianFromLatLongDeg)
  center <- rayProjectedMean(carts)
  geoLatLongDegFromCartesian(center)
}

#' Crude way to convert decimal latitudes and longitudes to approximate eastings and northings.
#' 
#' For serious applications you should probably use some combination of the rgdal and PBSmapping packages. Positive latitudes are north, while negative latitudes are south. Positive longitudes are east, while negative longitudes are west. Handles the prime antimeridian (near international date line) well. Does not work well in polar regions, where easting and northing become singular. This function produces approximate results, by assuming that the Earth is a perfect sphere of radius 6378 km. The output eastings and northings are relative to the directional mean of the locations, which you can using geoMeanLatLongDeg.
#' @param latLongDegs A list of 2D real vectors, where each vector is a latitude-longitude pair in decimal degrees.
#' @return A list of 2D real vectors. Easting-northing pairs, in km, relative to the mean location.
geoEastNorthsFromLatLongDegs <- function(latLongDegs) {
  carts <- lapply(latLongDegs, geoCartesianFromLatLongDeg)
  center <- rayProjectedMean(carts)
  up <- rayOrthogonalProjection(center, c(0, 0, 1))
  rot <- rbind(center, cross(up, center), up)
  tans <- lapply(carts, rayTangentVectorFromPoint, rot)
  lapply(tans, function(v) 6378 * v)
}

# Code to test those functions.
#latLongDegs <- list(c(0, 0), c(0, 15), c(10, 0))
#carts <- lapply(latLongDegs, geoCartesianFromLatLongDeg)
#plot3D(radius=NULL, points=carts)
#plot(x=sapply(latLongDegs, function(xy) xy[[2]]), y=sapply(latLongDegs, function(xy) xy[[1]]))
#eastNorths <- geoEastNorthsFromLatLongDegs(latLongDegs)
#plot(x=sapply(eastNorths, function(en) en[[1]]), y=sapply(eastNorths, function(en) en[[2]]))
#carts <- replicate(10, rayNormalized(mvrnorm(n=1, mu=c(-3, 0, 0), Sigma=diag(c(0.1, 0.1, 0.1)))), simplify=FALSE)
#plot3D(radius=NULL, points=carts)
#latLongDegs <- lapply(carts, geoLatLongDegFromCartesian)
#eastNorths <- geoEastNorthsFromLatLongDegs(latLongDegs)
#plot(x=sapply(eastNorths, function(en) en[[1]]), y=sapply(eastNorths, function(en) en[[2]]))



### PALEOMAGNETISM ###

# Given site lat-long, pmag trend-plunge (decl-incl), and pmag alpha95.
# Outputs paleopole lat-long, 95% confidence angles in lat-long.
# Notice that nothing is in degrees; everything is in radians.
# Based on appendix to Chapter 11 of Butler's 1992 textbook.
geoPaleoPoleFromDirection <- function(siteLat, siteLong, decl, incl, alpha95) {
  # Compute the paleopole latitude lambdap.
  p <- atan(2 / tan(incl))
  lambdap <- arcSin(sin(siteLat) * cos(p) + cos(siteLat) * sin(p) * cos(decl))
  # Compute the paleopole longitude phip.
  beta <- arcSin(sin(decl) * sin(p) / cos(lambdap))
  if (cos(p) >= sin(lambdap) * sin(siteLat))
    phip <- siteLong + beta
  else
    phip <- siteLong + pi - beta
  # Compute the 95% confidence angle for the paleopole latitude, in the usual crazy notation.
  dp <- 2 * alpha95 / (1 + 3 * cos(incl)^2)
  # Compute the 95% confidence angle for the paleopole longitude, in the usual crazy notation.
  dm <- alpha95 * sin(p) / cos(incl)
  # When the inclination is negative, Butler's algorithm seems to return paleo-south, not paleo-north.
  if (incl > 0)
    c(lambdap, phip %% (2 * pi), dp, dm)
  else
    c(-lambdap, (phip + pi) %% (2 * pi), dp, dm)
}

# Takes in A95, not alpha95.
geoPaleoDirectionFromPole <- function(siteLat, siteLong, poleLat, poleLong, a95) {
  p <- arcCos(sin(poleLat) * sin(siteLat) + cos(poleLat) * cos(siteLat) * cos(poleLong - siteLong))
  incl <- atan(2 / tan(p))
  dIncl <- 2 * a95 / (1 + 3 * cos(p)^2)
  decl <- arcCos((sin(poleLat) - sin(siteLat) * cos(p)) / (cos(siteLat) * sin(p)))
  # This correction is glossed over in Butler's Eq. (A.58).
  if ((siteLong - poleLong) %% (2 * pi) < pi)
    decl <- 2 * pi - decl
  dDecl <- arcSin(sin(a95) / sin(p))
  alpha95 <- arcSin(cos(incl) * sin(dDecl))
  c(decl, incl, dDecl, dIncl, alpha95)
}

# Test that they are inverses of each other.
geoTestPaleoPoleFromPole <- function() {
  siteLat <- runif(1, -pi / 2, pi / 2)
  siteLong <- runif(1, 0, 2 * pi)
  poleLat <- runif(1, -pi / 2, pi / 2)
  poleLong <- runif(1, 0, 2 * pi)
  dir <- geoPaleoDirectionFromPole(siteLat, siteLong, poleLat, poleLong, a95=(5 * degree))
  pole <- geoPaleoPoleFromDirection(siteLat, siteLong, dir[[1]], dir[[2]], alpha95=(5 * degree))
  c(siteLat, siteLong, dir[1:2], poleLat, poleLong, pole[1:2])
}
geoTestPaleoDirectionFromDirection <- function() {
  siteLat <- runif(1, -pi / 2, pi / 2)
  siteLong <- runif(1, 0, 2 * pi)
  decl <- runif(1, 0, 2 * pi)
  incl <- runif(1, -pi / 2, pi / 2)
  pole <- geoPaleoPoleFromDirection(siteLat, siteLong, decl, incl, alpha95=(5 * degree))
  dir <- geoPaleoDirectionFromPole(siteLat, siteLong, pole[[1]], pole[[2]], a95=(5 * degree))
  c(siteLat, siteLong, pole[1:2], c(decl, incl), dir[1:2])
}
#test <- replicate(1000, geoTestPaleoPoleFromPole())
#dot(test[5,] - test[7,], test[5,] - test[7,])
#dot(test[6,] - test[8,], test[6,] - test[8,])
#test <- replicate(1000, geoTestPaleoDirectionFromDirection())
#dot(test[5,] - test[7,], test[5,] - test[7,])
#dot(test[6,] - test[8,], test[6,] - test[8,])

# Other tests.
#geoPaleoPoleFromDirection(0 * degree, 0 * degree, 0 * degree, 45 * degree, 5 * degree) / degree
#geoPaleoPoleFromDirection(0 * degree, 0 * degree, 180 * degree, 45 * degree, 5 * degree) / degree
#geoPaleoPoleFromDirection(45 * degree, 0 * degree, 45 * degree, 60 * degree, 5 * degree) / degree



### DEDUCING ROTATION FROM PALEOMAGNETIC DIRECTION AND DIKE ASSUMPTION ###

#' Rotations deduced from two rays and a verticality assumption.
#' 
#' As in Allerton and Vine (1987). But to understand the mathematical treatment used here, see the appendix to our orientation statistics paper. We seek the rotation R that satisfies two constraints: R^T * dir = tmv and R^T * pole is horizontal. Based on these constraints, we get zero or two candidate rotations R at each station. In the zero case, we choose the unique rotation that comes closest to making R^T * pole horizontal.
#' @param pole A line (unit 3D real vector).
#' @param dir A ray (unit 3D real vector).
#' @param tmv A ray (unit 3D real vector).
#' @return A list of 3x3 real matrices (special orthogonal). The length of the list is 1 or 2.
geoRotationsFromDikePaleomag <- function(pole, dir, tmv) {
  # Build the rotations T, M from our appendix.
  v <- rayOrthogonalUniform(tmv)
  tMat <- rotProjectedMatrix(rbind(tmv, v, cross(tmv, v)))
  v <- rayOrthogonalUniform(dir)
  mMat <- rotProjectedMatrix(rbind(dir, v, cross(dir, v)))
  # Compute the coefficients a, b, c from our appendix.
  aa <- mMat[[1, 1]] * pole[[1]] * tMat[[1, 3]] + mMat[[1, 2]] * pole[[2]] * tMat[[1, 3]] +
    mMat[[1, 3]] * pole[[3]] * tMat[[1, 3]] - mMat[[2, 1]] * pole[[1]] * tMat[[2, 3]] -
    mMat[[2, 2]] * pole[[2]] * tMat[[2, 3]] - mMat[[2, 3]] * pole[[3]] * tMat[[2, 3]] -
    mMat[[3, 1]] * pole[[1]] * tMat[[3, 3]] - mMat[[3, 2]] * pole[[2]] * tMat[[3, 3]] -
    mMat[[3, 3]] * pole[[3]] * tMat[[3, 3]]
  bb <- 2 * mMat[[3, 1]] * pole[[1]] * tMat[[2, 3]] + 2 * mMat[[3, 2]] * pole[[2]] * tMat[[2, 3]] +
    2 * mMat[[3, 3]] * pole[[3]] * tMat[[2, 3]] - 2 * mMat[[2, 1]] * pole[[1]] * tMat[[3, 3]] -
    2 * mMat[[2, 2]] * pole[[2]] * tMat[[3, 3]] - 2 * mMat[[2, 3]] * pole[[3]] * tMat[[3, 3]]
  cc <- mMat[[1, 1]] * pole[[1]] * tMat[[1, 3]] + mMat[[1, 2]] * pole[[2]] * tMat[[1, 3]] +
    mMat[[1, 3]] * pole[[3]] * tMat[[1, 3]] + mMat[[2, 1]] * pole[[1]] * tMat[[2, 3]] +
    mMat[[2, 2]] * pole[[2]] * tMat[[2, 3]] + mMat[[2, 3]] * pole[[3]] * tMat[[2, 3]] +
    mMat[[3, 1]] * pole[[1]] * tMat[[3, 3]] + mMat[[3, 2]] * pole[[2]] * tMat[[3, 3]] +
    mMat[[3, 3]] * pole[[3]] * tMat[[3, 3]]
  # Solve for sigma and the final rotation.
  us <- realQuadraticSolutions(aa, bb, cc)
  if (length(us) >= 1) {
    sigmas <- sapply(us, function(u) atan2(2 * u, 1 - u^2))
    rs <- lapply(sigmas, function(sigma) {t(mMat) %*% rotMatrixAboutX(sigma) %*% tMat})
  } else {
    alpha <- mMat[[2, 1]] * pole[[1]] * tMat[[2, 3]] + mMat[[2, 2]] * pole[[2]] * tMat[[2, 3]] +
      mMat[[2, 3]] * pole[[3]] * tMat[[2, 3]] + mMat[[3, 1]] * pole[[1]] * tMat[[3, 3]] +
      mMat[[3, 2]] * pole[[2]] * tMat[[3, 3]] + mMat[[3, 3]] * pole[[3]] * tMat[[3, 3]]
    beta <- mMat[[3, 1]] * pole[[1]] * tMat[[2, 3]] + mMat[[3, 2]] * pole[[2]] * tMat[[2, 3]] +
      mMat[[3, 3]] * pole[[3]] * tMat[[2, 3]] - mMat[[2, 1]] * pole[[1]] * tMat[[3, 3]] -
      mMat[[2, 2]] * pole[[2]] * tMat[[3, 3]] - mMat[[2, 3]] * pole[[3]] * tMat[[3, 3]]
    sigmas <- c(atan2(beta, alpha), atan2(beta, alpha) + pi)
    rs <- lapply(sigmas, function(sigma) {t(mMat) %*% rotMatrixAboutX(sigma) %*% tMat})
    # On 2019/07/19 Sarah made me realize that r should be t(r) in the next line.
    absZs <- lapply(rs, function(r) abs((t(r) %*% pole)[[3]]))
    rs <- rs[which.min(absZs)]
  }
  rs
}



### LOADING DATA FROM FILE ###

# Helper function for geoDataFromFile.
geoPoleDirectionRotationFromData <- function(dataFrame) {
  newFrame <- emptyDataFrame(nrow(dataFrame))
  if (!is.null(dataFrame$strikeDeg) && !is.null(dataFrame$dipDeg))
    if (!is.null(dataFrame$rakeDeg)) {
      # Make rotation, pole, direction from strike, dip, rake.
      newFrame$rotation <- lapply(1:nrow(dataFrame),
                                  function(i) geoRotationFromStrikeDipRakeDeg(
                                    c(dataFrame[i,]$strikeDeg, dataFrame[i,]$dipDeg, dataFrame[i,]$rakeDeg)))
      newFrame$pole <- lapply(1:nrow(newFrame), function(i) newFrame$rotation[[i]][1,])
      newFrame$direction <- lapply(1:nrow(newFrame), function(i) newFrame$rotation[[i]][2,])
    } else if (!is.null(dataFrame$trendDeg) && !is.null(dataFrame$plungeDeg)) {
      # Make rotation, pole, direction from strike, dip, trend, plunge.
      newFrame$pole <- thread(
        function(strikeDeg, dipDeg) geoCartesianFromStrikeDipDeg(c(strikeDeg, dipDeg)),
        dataFrame$strikeDeg, dataFrame$dipDeg)
      newFrame$direction <- thread(
        function(trendDeg, plungeDeg) geoCartesianFromTrendPlungeDeg(c(trendDeg, plungeDeg)),
        dataFrame$trendDeg, dataFrame$plungeDeg)
      newFrame$rotation <- thread(
        function(pole, direction) rotProjectedMatrix(rbind(pole, direction, cross(pole, direction))),
        newFrame$pole, newFrame$direction)
    } else {
      # Make pole from strike, dip.
      newFrame$pole <- lapply(1:nrow(dataFrame),
                              function(i) geoCartesianFromStrikeDipDeg(c(dataFrame[i,]$strikeDeg, dataFrame[i,]$dipDeg)))
    }
  else if (!is.null(dataFrame$trendDeg) && !is.null(dataFrame$plungeDeg)) {
    # Make direction from trend, plunge.
    newFrame$direction <- lapply(1:nrow(dataFrame),
                                 function(i) geoCartesianFromTrendPlungeDeg(c(dataFrame[i,]$trendDeg, dataFrame[i,]$plungeDeg)))
  }
  if (!is.null(dataFrame$trend1Deg) && !is.null(dataFrame$plunge1Deg) &&
        !is.null(dataFrame$trend2Deg) && !is.null(dataFrame$plunge2Deg)) {
    # Make rotation from two trends and plunges.
    f <- function(i) {
      axis1 <- geoCartesianFromTrendPlungeDeg(c(dataFrame[i,]$trend1Deg, dataFrame[i,]$plunge1Deg))
      axis2 <- geoCartesianFromTrendPlungeDeg(c(dataFrame[i,]$trend2Deg, dataFrame[i,]$plunge2Deg))
      rotProjectedMatrix(rbind(axis1, axis2, cross(axis1, axis2)))
    }
    newFrame$rotation <- lapply(1:nrow(dataFrame), f)
  }
  if (!is.null(dataFrame$x1) && !is.null(dataFrame$y1) && !is.null(dataFrame$z1) &&
        !is.null(dataFrame$x2) && !is.null(dataFrame$y2) && !is.null(dataFrame$z2)) {
    # Make rotation from its first two rows.
    f <- function(i) {
      axis1 <- c(dataFrame[i,]$x1, dataFrame[i,]$y1, dataFrame[i,]$z1)
      axis2 <- c(dataFrame[i,]$x2, dataFrame[i,]$y2, dataFrame[i,]$z2)
      rotProjectedMatrix(rbind(axis1, axis2, cross(axis1, axis2)))
    }
    newFrame$rotation <- lapply(1:nrow(dataFrame), f)
  }
  if (!is.null(dataFrame$x) && !is.null(dataFrame$y) && !is.null(dataFrame$z)) {
    # Make direction from x, y, z Cartesian coordinates.
    newFrame$direction <- lapply(1:nrow(dataFrame), function(i) unit(c(dataFrame[i,]$x, dataFrame[i,]$y, dataFrame[i,]$z)))
  }
  newFrame
}

# Helper function for geoDataFromFile, etc.
# Assumes that a1, a2, a3, rotation fields present in dataFrame. Already normalized if they are supposed to be.
# Returns a new data frame with fields a, logA, tensor, vector.
geoEllipsoidFromRotationAData <- function(dataFrame, isNormalized) {
  newFrame <- emptyDataFrame(nrow(dataFrame))
  # Bind semi-axis lengths a1, a2, a3 into arrays.
  newFrame$a <- lapply(1:nrow(dataFrame), function(i) c(dataFrame$a1[[i]], dataFrame$a2[[i]], dataFrame$a3[[i]]))
  newFrame$logA <- lapply(1:nrow(dataFrame), function(i) log(newFrame$a[[i]]))
  # Make ellipsoid tensor from a1, a2, a3, rotation.
  newFrame$tensor <- lapply(1:nrow(dataFrame),
                            function(i) {t(dataFrame$rotation[[i]]) %*%
                                           diag(c(dataFrame$a1[[i]], dataFrame$a2[[i]], dataFrame$a3[[i]])^-2) %*%
                                           dataFrame$rotation[[i]]})
  # Make log-ellipsoid tensor, but don't bother storing it.
  logElls <- lapply(1:nrow(dataFrame),
                    function(i) {t(dataFrame$rotation[[i]]) %*% diag(-2 * newFrame$logA[[i]]) %*% dataFrame$rotation[[i]]})
  # Make ellipsoid vector from log-ellipsoid tensor.
  f <- function(i) {
    if (isNormalized)
      ellNormalizedVectorFromLog(logElls[[i]])
    else
      ellVectorFromLog(logElls[[i]])
  }
  newFrame$vector <- lapply(1:nrow(newFrame), f)
  newFrame
}

# Helper function for geoDataFromFile.
geoEllipsoidFromVectorData <- function(dataFrame) {
  newFrame <- emptyDataFrame(nrow(dataFrame))
  f <- function(i) {
    v <- c(dataFrame$vector1[[i]], dataFrame$vector2[[i]], dataFrame$vector3[[i]], dataFrame$vector4[[i]], dataFrame$vector5[[i]])
    if (is.null(dataFrame$vector6))
      v <- c(v, dataFrame$vector6[[i]])
    ell <- ellEllipsoidFromVector(v)
    }
  ellipsoids <- lapply(1:nrow(dataFrame), f)
  newFrame$vector <- lapply(ellipsoids, function(ell) ell$vector)
  newFrame$rotation <- lapply(ellipsoids, function(ell) ell$rotation)
  newFrame$tensor <- lapply(ellipsoids, function(ell) ell$tensor)
  newFrame$a <- lapply(ellipsoids, function(ell) ell$a)
  newFrame$logA <- lapply(ellipsoids, function(ell) ell$logA)
  newFrame
}

#' Loads the file as a data frame. Attempts to infer pole, direction, rotation, and/or ellipsoid from each datum.
#'
#' In addition to loading the file into an R data frame (table), this function appends "Deg" to the names of certain strike, dip, trend, and plunge fields. It also attempts to infer poles, directions, rotations, and ellipsoids from the data, according to the following rules. When two rules conflict, the behavior is undefined. So avoid conflicts yourself.
#' If fields strike-dip exist, then adds a field pole. If field rake also exists, then adds a field direction.
#' If fields trend-plunge or x-y-z exist, then adds a field direction.
#' If fields trend1-plunge1-trend2-plunge2 or pole-direction or x1-y1-z1-x2-y2-z2, then adds a field rotation.
#' If fields a1-a2, but not field a3, then adds a field a3 such that a1 a2 a3 == 1.
#' If fields rotation-a1-a2-a3, then adds ellipsoid fields a, logA, tensor, vector.
#' If fields vector1-vector2-vector3-vector4-vector5(-vector6), then adds ellipsoid fields a, logA, rotation, tensor, vector.
#' vector is 5D if normalized, 6D if not.
#' @param fileName Character. The name of the file (or path relative to the working directory).
#' @param doNormalize Logical. Whether to normalize ellipsoids to have volume equal to that of the unit sphere.
#' @param separator Character or NULL. The separator for fields in the file. If NULL, then attempts to deduce either ',' or '\t' from the file name suffix.
#' @return Data frame. All fields from the original file, maybe some with 'Deg' appended, and maybe a bunch of other fields as described above.
geoDataFromFile <- function(fileName, doNormalize=FALSE, separator=NULL) {
  # Guess the file format, defaulting to tab-separated values.
  if (is.null(separator)) {
    if (substr(fileName, nchar(fileName) - 3, nchar(fileName)) == ".csv")
      separator <- ","
    else
      separator <- "\t"
  }
  # Read the file into a data frame.
  dataFrame <- read.table(fileName, header=TRUE, sep=separator)
  # Append "Deg" to certain column names.
  for (name in c("strike", "dip", "rake", "trend", "plunge", "trend1", "plunge1", "trend2", "plunge2"))
    if (length(which(names(dataFrame) == name)) == 1)
      names(dataFrame)[[which(names(dataFrame) == name)]] <- paste0(name, "Deg")
  # Post-process pole, direction, rotation.
  dataFrame <- cbind(dataFrame, geoPoleDirectionRotationFromData(dataFrame))
  # Post-process ellipsoid semi-axes based on a1 a2 a3 == 1.
  if (!is.null(dataFrame$a1) && !is.null(dataFrame$a2)) {
    if (is.null(dataFrame$a3))
      dataFrame$a3 <- sapply(1:nrow(dataFrame), function(i) {1 / (dataFrame$a1[[i]] * dataFrame$a2[[i]])})
    else if (doNormalize) {
      scalars <- sapply(1:nrow(dataFrame), function(i) (dataFrame$a1[[i]] * dataFrame$a2[[i]] * dataFrame$a3[[i]])^(-1 / 3))
      dataFrame$a1 <- dataFrame$a1 * scalars
      dataFrame$a2 <- dataFrame$a2 * scalars
      dataFrame$a3 <- dataFrame$a3 * scalars
    }
  }
  # Post-process ellipsoid and its auxiliary information.
  if (!is.null(dataFrame$a1) && !is.null(dataFrame$a2) && !is.null(dataFrame$a3) && !is.null(dataFrame$rotation))
    dataFrame <- cbind(dataFrame, geoEllipsoidFromRotationAData(dataFrame, doNormalize))
  if (!is.null(dataFrame$vector1) && !is.null(dataFrame$vector2) && !is.null(dataFrame$vector3) && 
        !is.null(dataFrame$vector4) && !is.null(dataFrame$vector5))
    dataFrame <- cbind(dataFrame, geoEllipsoidFromVectorData(dataFrame))
  dataFrame
}



### SPECIALIZED READING DATA FROM FILE ###

# Helper function for geoDataFromAGICOFile.
geoSemiaxisDataFromAGICOData <- function(dataFrame, meanSuscepts, doNormalize) {
  f <- function(i) {
    k <- rbind(c(dataFrame[i,]$K11, dataFrame[i,]$K12, dataFrame[i,]$K13),
          c(dataFrame[i,]$K12, dataFrame[i,]$K22, dataFrame[i,]$K23),
          c(dataFrame[i,]$K13, dataFrame[i,]$K23, dataFrame[i,]$K33))
    eig <- eigen(k, symmetric=TRUE, only.values=TRUE)$values
    eig <- eig * meanSuscepts[[i]]
    if (doNormalize)
      eig <- eig / prod(eig)^(1 / 3)
    eig
  }
  eigs <- lapply(1:nrow(dataFrame), f)
  data.frame(
    a1=sapply(eigs, function(eig) eig[[1]]),
    a2=sapply(eigs, function(eig) eig[[2]]),
    a3=sapply(eigs, function(eig) eig[[3]]))
}

#' Loading a file of AMS ellipsoids in a certain AGICO file format.
#' 
#' Warning: doNormalize=FALSE doesn't work?!! In addition to an AGICO-processed kappa bridge data file listing AMS ellipsoids, you need the mean susceptibility for each AMS ellipsoid, which is output in a separate file. The mean susceptibilities tell us how the ellipsoids were normalized, so that we can de-normalize them. Warning: The conventions used by various authors for AMS ellipsoids are so confusing and contradictory that I'm not completely sure that the ellipsoid semi-axis lengths are configured correctly. The colloquial term for this kind of situation is 'Hell on Earth'.
#' @param fileName Character. The name of the file (or path relative to the working directory).
#' @param meanSuscepts A vector of real numbers. The mean susceptibilities corresponding to the rows in the data file.
#' @param doNormalize Logical. Whether to normalize ellipsoids to have volume equal to that of the unit sphere.
#' @param separator Character or NULL. The separator for fields in the file. If NULL, then attempts to deduce either ',' or '\t' from the file name suffix.
#' @return Data frame. All fields from the original file, maybe with a bunch of ellipsoid fields appended.
geoEllipsoidDataFromAGICOFile <- function(fileName, meanSuscepts, doNormalize=FALSE, separator=NULL) {
  # Guess the file format, defaulting to tab-separated values.
  if (is.null(separator)) {
    if (substr(fileName, nchar(fileName) - 3, nchar(fileName)) == ".csv")
      separator <- ","
    else
      separator <- "\t"
  }
  # Read the file into a data frame.
  dataFrame <- read.table(fileName, header=TRUE, sep=separator)
  # Infer the ellipsoid semi-axis lengths.
  dataFrame <- cbind(
    dataFrame, 
    geoSemiaxisDataFromAGICOData(dataFrame, meanSuscepts, doNormalize))
  # Build the rotation field.
  p1s <- lapply(
    1:nrow(dataFrame), 
    function(i) geoCartesianFromTrendPlungeDeg(c(dataFrame$K1dec[[i]], dataFrame$K1inc[[i]])))
  p2s <- lapply(
    1:nrow(dataFrame), 
    function(i) geoCartesianFromTrendPlungeDeg(c(dataFrame$K2dec[[i]], dataFrame$K2inc[[i]])))
  dataFrame$rotation <- lapply(
    1:nrow(dataFrame), 
    function(i) rotProjectedMatrix(rbind(p1s[[i]], p2s[[i]], cross(p1s[[i]], p2s[[i]]))))
  # Finish building the ellipsoid fields from a1-a2-a3-rotation.
  dataFrame <- cbind(dataFrame, geoEllipsoidFromRotationAData(dataFrame, doNormalize))
  dataFrame
}

#' Loading a file of AMS ellipsoids in a certain Institute for Rock Magnetism (IRM) file format.
#' 
#' Dmax, Imax, and max are the trend (declination), plunge (inclination), and length of the longest semi-axis of the magnitude ellipsoid. Length is also the largest eigenvalue of the susceptibility tensor (see e.g. Hrouda, 1982). Similarly, Dmin, Imin, min are for the shortest semi-axis. (This is all according to Mike Jackson, personal communication, IRM, 2016/03/09. He told me that the max, int, min in the file are the eigenvalues of the susceptibility tensor, which are also the semi-axis lengths of the magnitude ellipsoid (Hrouda, 1982), which is the ellipsoid that all recent authors care about.)
#' @param fileName Character. The name of the file (or path relative to the working directory).
#' @param doNormalize Logical. Whether to normalize ellipsoids to have volume equal to that of the unit sphere.
#' @param separator Character or NULL. The separator for fields in the file. If NULL, then attempts to deduce either ',' or '\t' from the file name suffix.
#' @return Data frame. All fields from the original file, maybe with a bunch of ellipsoid fields appended.
geoEllipsoidDataFromIRMFile <- function(fileName, doNormalize=FALSE, separator=NULL) {
  # Guess the file format, defaulting to tab-separated values.
  if (is.null(separator)) {
    if (substr(fileName, nchar(fileName) - 3, nchar(fileName)) == ".csv")
      separator <- ","
    else
      separator <- "\t"
  }
  # Read the file into a data frame.
  dataFrame <- read.table(fileName, header=TRUE, sep=separator)
  # Make a new data frame with rotation field.
  newFrame <- emptyDataFrame(nrow(dataFrame))
  newFrame$trend1Deg <- dataFrame$Dmin..in.situ.
  newFrame$plunge1Deg <- dataFrame$Imin..in.situ.
  newFrame$trend2Deg <- dataFrame$Dmax..in.situ.
  newFrame$plunge2Deg <- dataFrame$Imax..in.situ.
  f <- function(i) {
    axis1 <- geoCartesianFromTrendPlungeDeg(c(newFrame$trend1Deg[[i]], newFrame$plunge1Deg[[i]]))
    axis2 <- geoCartesianFromTrendPlungeDeg(c(newFrame$trend2Deg[[i]], newFrame$plunge2Deg[[i]]))
    rotProjectedMatrix(rbind(axis1, axis2, cross(axis1, axis2)))
  }
  newFrame$rotation <- lapply(1:nrow(newFrame), f)
  # Add fields a1, a2, a3, possibly normalized.
  newFrame$a1 <- dataFrame$min
  newFrame$a2 <- dataFrame$max
  newFrame$a3 <- dataFrame$int
  if (doNormalize) {
    f <- function(i) {
      (newFrame$a1[[i]] * newFrame$a2[[i]] * newFrame$a3[[i]])^(-1 / 3)
    }
    scalars <- sapply(1:nrow(newFrame), f)
    newFrame$a1 <- newFrame$a1 * scalars
    newFrame$a2 <- newFrame$a2 * scalars
    newFrame$a3 <- newFrame$a3 * scalars
  }
  # Add ellipsoid fields: a, logA, tensor, vector.
  newFrame <- cbind(newFrame, geoEllipsoidFromRotationAData(newFrame, doNormalize))
  dataFrame <- cbind(dataFrame, newFrame)
  dataFrame
}

#' Loading ellipsoids from a certain Avizo file format.
#' 
#' The Avizo image analysis software is used to process some X-ray computed tomography data. In 2015 Vasilis Chatzaras and I were confused about whether the EigenVal fields were semi-axis lengths or some other power of them. Avizo tech support stated that they were in units of length squared. That is, they are roughly proportional to semi-axis lengths squared. But I think what's really going on is this: The software records the (x, y, z) coordinates of each voxel. Then it computes the 3x3 covariance matrix of that {(x, y, z)_1, ..., (x, y, z)_n} data set. Then it reports the eigensystem of that matrix. Testing on 2017 August 20 by Vasilis, with data sets where the actual extent of the objects is known, suggests that the actual extent is 3.5 to 4.5 times as long as the square roots of the eigenvalues. That's a large amount of slop. Hmm.
#' @param fileName Character. The name of the file (or path relative to the working directory).
#' @param doNormalize Logical. Whether to normalize ellipsoids to have volume equal to that of the unit sphere.
#' @param separator Character or NULL. The separator for fields in the file. If NULL, then attempts to deduce either ',' or '\t' from the file name suffix.
#' @return Data frame. All fields from the original file, maybe with a bunch of ellipsoid fields appended.
geoEllipsoidDataFromAvizoFile <- function(fileName, doNormalize=FALSE, separator=NULL) {
  # Guess the file format, defaulting to tab-separated values.
  if (is.null(separator)) {
    if (substr(fileName, nchar(fileName) - 3, nchar(fileName)) == ".csv")
      separator <- ","
    else
      separator <- "\t"
  }
  # Read the file into a data frame.
  dataFrame <- read.table(fileName, header=TRUE, sep=separator)
  # Rename the eigenvector columns.
  names(dataFrame)[[which(names(dataFrame) == "EigenVec1X")]] <- "x1"
  names(dataFrame)[[which(names(dataFrame) == "EigenVec1Y")]] <- "y1"
  names(dataFrame)[[which(names(dataFrame) == "EigenVec1Z")]] <- "z1"
  names(dataFrame)[[which(names(dataFrame) == "EigenVec2X")]] <- "x2"
  names(dataFrame)[[which(names(dataFrame) == "EigenVec2Y")]] <- "y2"
  names(dataFrame)[[which(names(dataFrame) == "EigenVec2Z")]] <- "z2"
  # The eigenvalues are the semi-axis lengths squared.
  dataFrame$a1 <- sqrt(dataFrame$EigenVal1)
  dataFrame$a2 <- sqrt(dataFrame$EigenVal2)
  dataFrame$a3 <- sqrt(dataFrame$EigenVal3)
  # Add the rotation field.
  dataFrame <- cbind(dataFrame, geoPoleDirectionRotationFromData(dataFrame))
  # Add the ellipsoid vector, tensor, a, logA fields.
  if (doNormalize) {
    scalars <- sapply(1:nrow(dataFrame), function(i) (dataFrame$a1[[i]] * dataFrame$a2[[i]] * dataFrame$a3[[i]])^(-1 / 3))
    dataFrame$a1 <- dataFrame$a1 * scalars
    dataFrame$a2 <- dataFrame$a2 * scalars
    dataFrame$a3 <- dataFrame$a3 * scalars
  }
  # Post-process ellipsoid and its auxiliary information.
  dataFrame <- cbind(dataFrame, geoEllipsoidFromRotationAData(dataFrame, doNormalize))
  dataFrame
}



### WRITING DATA TO FILE ###

# Not well tested. Because R's lists are not really first-class objects? Because R's basic data structures are insanely designed?
geoFileFromData <- function(fileName, dataFrame, separator=NULL) {
  # Guess the file format, defaulting to tab-separated values.
  if (is.null(separator)) {
    if (substr(fileName, nchar(fileName) - 3, nchar(fileName)) == ".csv")
      separator <- ","
    else
      separator <- "\t"
  }
  # Write it.
  write.table(dataFrame, file=fileName, sep=separator, row.names=FALSE, quote=FALSE)
}

# Helper function for geoFaultKinFileFromMatrices.
geoFaultKinFromPoleVorticity <- function(r) {
  q <- geoPoleHangingFromPoleVorticity(r)
  sdDeg <- geoStrikeDipDegFromCartesian(q[1,])
  tpDeg <- geoTrendPlungeDegFromCartesian(q[2,])
  if (tpDeg[[2]] >= 0)
    sense <- 999999
  else {
    sense <- -999999
    tpDeg <- c((tpDeg[[1]] + 180) %% 360, -tpDeg[[2]])
  }
  c(sdDeg, tpDeg, sense)
}

#' Writing a file of faults-with-slip, in a certain FaultKin format.
#' 
#' Writes a tab-separated values (TSV) file to disk. You should use the suffix '.txt', not '.tsv'. Then you can pretty easily load this file into the FaultKin software by Allmendinger and Cardozo.
#' @param rs A list of rotation matrices. The orientations of the faults-with-slip. Each matrix has the fault pole as its first row, the vorticity vector of slip as its second row, and the cross product of those as its third row.
#' @param fileName Character. The name of the file (or path relative to the working directory). It should end in the suffix '.txt'.
#' @return NULL.
geoFaultKinFileFromMatrices <- function(rs, fileName) {
  mat <- t(sapply(rs, geoFaultKinFromPoleVorticity))
  df <- emptyDataFrame(nrow(mat))
  df[["Fault strike"]] <- mat[,1]
  df[["Fault dip"]] <- mat[,2]
  df[["Straie trend"]] <- mat[,3]
  df[["Straie dip"]] <- mat[,4]
  df[["Sense of slip"]] <- sapply(mat[,5], function(s) {if (s < 0) "T" else "N"})
  write.table(df, file=fileName, sep="\t", row.names=FALSE, quote=FALSE,
              col.names=c("Fault strike", "Fault dip", "Straie trend", "Straie plunge", "Sense of Slip"))
}

#' Writing a file of trends and plunges in degrees.
#' 
#' Writes a tab-separated values (TSV) file to disk. You can use whatever suffix you like, but the obvious choices are '.tsv' and '.txt'. The file has two columns, headered 'trend' and 'plunge'.
#' @param us A list of rays (length-1 3D vectors).
#' @param fileName Character. The name of the file (or path relative to the working directory).
#' @return NULL.
geoTrendPlungeDegFileFromRays <- function(us, fileName) {
  tpDegs <- lapply(us, geoTrendPlungeDegFromCartesian)
  df <- data.frame(trend=sapply(tpDegs, function(tpDeg) tpDeg[[1]]),
                   plunge=sapply(tpDegs, function(tpDeg) tpDeg[[2]]))
  write.table(df, file=fileName, sep="\t", row.names=FALSE, quote=FALSE)
}



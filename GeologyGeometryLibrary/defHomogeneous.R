


# Copyright 2016-2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### HOMOGENEOUS SIMPLE SHEAR ###

#' Velocity gradient tensor for homogeneous simple shear along the x-z-plane.
#' 
#' @param gamma A real number. Positive for dextral shear, negative for sinistral shear.
#' @return A real 3x3 matrix.
defSimpleVGT <- function(gamma) {
  rbind(c(0, gamma, 0), c(0, 0, 0), c(0, 0, 0))
}

#' Position gradient tensor for homogeneous simple shear along the x-z-plane.
#' 
#' @param gamma A real number. Positive for dextral shear, negative for sinistral shear.
#' @return A real 3x3 matrix.
defSimplePGT <- function(gamma) {
  rbind(c(1, gamma, 0), c(0, 1, 0), c(0, 0, 1))
}



### HOMOGENEOUS MONOCLINIC TRANSPRESSION (FOSSEN AND TIKOFF, 1993) ###

#' Velocity gradient tensor for homogeneous monoclinic transpression along the x-z-plane.
#' 
#' See Fossen and Tikoff (1993).
#' @param gamma A real number. Positive for dextral shear, negative for sinistral shear.
#' @param logK A real number. Negative for transpression, positive for transtension.
#' @return A real 3x3 matrix.
defMonoclinicVGT <- function(gamma, logK) {
  matrix(c(0, 0, 0, gamma, logK, 0, 0, 0, -logK), 3, 3)
}

#' Position gradient tensor for homogeneous monoclinic transpression along the x-z-plane.
#' 
#' See Fossen and Tikoff (1993).
#' @param gamma A real number. Positive for dextral shear, negative for sinistral shear.
#' @param logK A real number. Negative for transpression, positive for transtension.
#' @return A real 3x3 matrix.
defMonoclinicPGT <- function(gamma, logK) {
  if (logK == 0)
    rbind(c(1, gamma, 0), c(0, 1, 0), c(0, 0, 1))
  else {
    k <- exp(logK)
    rbind(c(1, gamma * (k - 1) / logK, 0), c(0, k, 0), c(0, 0, 1 / k))
  }
}

#' Inverse position gradient tensor for homogeneous monoclinic transpression along the x-z-plane.
#' 
#' The output is the inverse of the output of defMonoclinicPGT.
#' @param gamma A real number. Positive for dextral shear, negative for sinistral shear.
#' @param logK A real number. Negative for transpression, positive for transtension.
#' @return A real 3x3 matrix.
defMonoclinicInversePGT <- function(gamma, logK) {
  if (logK == 0)
    rbind(c(1, -gamma, 0), c(0, 1, 0), c(0, 0, 1))
  else {
    k <- exp(logK)
    rbind(c(1, gamma * (1 - k) / (k * logK), 0), c(0, 1 / k, 0), c(0, 0, k))
  }
}

#' Critical curve for homogeneous monoclinic transpression.
#' 
#' For certain combinations of (gamma, logK), homogeneous monoclinic transpression produces a finite strain ellipsoid that is an oblate spheroid. Thus the long axis direction, which is frequently interpreted as the lineation direction, is undefined. For any logK, this function returns the corresponding gamma >= 0 to make this happen.
#' @param logK A real number. Negative for transpression, positive for transtension.
#' @return A real number, non-negative. The gamma corresponding to logK.
defMonoclinicCriticalGammaFromLogK <- function(logK) {
  k <- exp(logK)
  -(k + 1) * sqrt(k^2 + 1) * logK / k
}

#' Angle of oblique convergence for homogeneous monoclinic transpression.
#' 
#' This is the horizontal angle between the movement direction and the shear zone. For dextral transpression, it is between 0 and pi / 2. Together with magnitude, it forms polar coordinates on the space of homogeneous monoclinic transpressions.
#' @param gamma A real number. Positive for dextral shear, negative for sinistral shear.
#' @param logK A real number. Negative for transpression, positive for transtension.
#' @return A real number. The angle of oblique convergence, in radians.
defMonoclinicAOCFromGammaLogK <- function(gamma, logK) {
  atan2(-logK, gamma)
}

#' Magnitude of homogeneous monoclinic transpression.
#' 
#' Because of the nature of steady homogeneous deformations, the magnitude can be interpreted in various ways. On the one hand, it measures the intensity of a deformation run for a fixed time interval, say from t = 0 to t = 1. On the other hand, it measures the amount of time that a deformation runs, if we fix that deformation's intensity. Anyway, this is a good measure of 'how far away from no deformation at all' the transpression is. Together with angle of oblique convergence, it forms polar coordinates on the space of homogeneous monoclinic transpressions.
#' @param gamma A real number. Positive for dextral shear, negative for sinistral shear.
#' @param logK A real number. Negative for transpression, positive for transtension.
#' @return A real number. The angle of oblique convergence, in radians.
defMonoclinicMagnitudeFromGammaLogK <- function(gamma, logK) {
  sqrt(gamma^2 + logK^2)
}

#' Parameter gamma for homogeneous monoclinic transpression, from other parameters.
#' 
#' This function and defMonoclinicLogKFromAOCMagnitude are inverse to defMonoclinicCriticalGammaFromLogK and defMonoclinicAOCFromGammaLogK.
#' @param aoc A real number.  The angle of oblique convergence, in radians.
#' @param mag A real number. The magnitude.
#' @return A real number, gamma.
defMonoclinicGammaFromAOCMagnitude <- function(aoc, mag) {
  mag * cos(aoc)
}

#' Parameter log k for homogeneous monoclinic transpression, from other parameters.
#' 
#' This function and defMonoclinicGammaFromAOCMagnitude are inverse to defMonoclinicCriticalGammaFromLogK and defMonoclinicAOCFromGammaLogK.
#' @param aoc A real number.  The angle of oblique convergence, in radians.
#' @param mag A real number. The magnitude.
#' @return A real number, logK.
defMonoclinicLogKFromAOCMagnitude <- function(aoc, mag) {
  -mag * sin(aoc)
}



### HOMOGENEOUS TRICLINIC TRANSPRESSION (JONES AND HOLDSWORTH, 1998; LIN ET AL., 1998) ###

#' Velocity gradient tensor for homogeneous triclinic transpression along a vertical, EW-striking shear plane.
#' 
#' This deformation is homogeneous, so there are effectively only two positions: inside the shear zone versus outside. If you just pass the vector u, then the shear zone is assumed to be infinitely thick, so the station in question is inside the shear zone. If you pass any more parameters than u, then pass all of them. They are used to determine whether the station is inside the shear zone at time t or not.
#' @param u A 3-dimensional real vector. The movement vector of the bounding rigid blocks relative to the shear plane.
#' @param hOf1 A real number, positive. The final half-width of the shear zone.
#' @param xOfT A 3-dimensional real vector. The location of the point, at which we wish to compute the VGT. Because the deformation is homogeneous, all points inside the zone have the same VGT, and all points outside have VGT = 0.
#' @param tt A real number. The time t at which x(t) = xOfT.
#' @return A real 3x3 matrix.
defTriclinicVGT <- function(u, hOf1=1, xOfT=c(0, 0, 0), tt=1) {
  hOfT <- hOf1 * exp(0.5 * u[[2]] * (tt - 1));
  if (xOfT[[2]]^2 >= hOfT^2)
    diag(c(0, 0, 0))
  else
    cbind(c(0, 0, 0), u, c(0, 0, -u[[2]]))
}

#' Position gradient tensor for homogeneous triclinic transpression along a vertical, EW-striking shear plane.
#' 
#' Equals defExp(defTriclinicVGT(...)), but faster and more robustly. Same parameters as defTriclinicVGT.!!
defTriclinicPGT <- function(u, hOf1=1, xOfT=c(0, 0, 0), tt=1) {
  hOfT <- hOf1 * exp(0.5 * u[[2]] * (tt - 1));
  if (xOfT[[2]]^2 >= hOfT^2)
    diag(c(1, 1, 1))
  else if (u[[2]] == 0)
    cbind(c(1, 0, 0), c(u[[1]], 1, u[[3]]), c(0, 0, 1))
  else {
    expU2 <- exp(u[[2]])
    rbind(
      c(1, (expU2 - 1) * u[[1]] / u[[2]], 0),
      c(0, expU2, 0),
      c(0, (expU2 - 1 / expU2) * u[[3]] / (2 * u[[2]]), 1 / expU2))
  }
}



### HOMOGENEOUS TRANSPORT TRANSPRESSION ###

#' Velocity gradient tensor for homogeneous 'transport' transpression along a vertical, EW-striking shear plane.
#' 
#' This deformation is homogeneous, so there are effectively only two positions: inside the shear zone versus outside. If you just pass the vector u, then the shear zone is assumed to be infinitely thick, so the station in question is inside the shear zone. If you pass any more parameters than u, then pass all of them. They are used to determine whether the station is inside the shear zone at time t or not.
#' @param u A 3-dimensional real vector. The movement vector of the bounding rigid blocks relative to the shear plane.
#' @param hOf1 A real number, positive. The final half-width of the shear zone.
#' @param xOfT A 3-dimensional real vector. The location of the point, at which we wish to compute the VGT. Because the deformation is homogeneous, all points inside the zone have the same VGT, and all points outside have VGT = 0.
#' @param tt A real number. The time t at which x(t) = xOfT.
#' @return A real 3x3 matrix.
defTransportVGT <- function(u, hOf1=1, xOfT=c(0, 0, 0), tt=1) {
  hOfT <- hOf1 * exp(0.5 * u[[2]] * (tt - 1));
  if (xOfT[[2]]^2 >= hOfT^2)
    diag(c(0, 0, 0))
  else
    rbind(
      c(0, u[[1]], 0),
      c(0, u[[2]], 0),
      c(0, 2 * u[[3]], -u[[2]]))
}

#' Position gradient tensor for homogeneous 'transport' transpression along a vertical, EW-striking shear plane.
#' 
#' Equals defExp(defTransportVGT(...)), but more quickly and robustly. The parameters are the same as in defTransportVGT.
#' @param u A 3-dimensional real vector. The movement vector of the bounding rigid blocks relative to the shear plane.
#' @param hOf1 A real number, positive. The final half-width of the shear zone.
#' @param xOfT A 3-dimensional real vector. The location of the point, at which we wish to compute the VGT. Because the deformation is homogeneous, all points inside the zone have the same VGT, and all points outside have VGT = 0.
#' @param tt A real number. The time t at which x(t) = xOfT.
#' @return A real 3x3 matrix.
defTransportPGT <- function(u, hOf1=1, xOfT=c(0, 0, 0), tt=1) {
  hOfT <- hOf1 * exp(0.5 * u[[2]] * (tt - 1));
  if (xOfT[[2]]^2 >= hOfT^2)
    diag(c(1, 1, 1))
  else if (u[[2]] == 0)
    cbind(c(1, 0, 0), c(u[[1]], 1, 2 * u[[3]]), c(0, 0, 1))
  else {
    expU2 <- exp(u[[2]])
    rbind(
      c(1, (expU2 - 1) * u[[1]] / u[[2]], 0),
      c(0, expU2, 0),
      c(0, (expU2 - 1 / expU2) * u[[3]] / u[[2]], 1 / expU2))
  }
}



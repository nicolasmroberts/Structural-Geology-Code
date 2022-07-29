


# Copyright 2016-2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This file offers a few functions dealing with deformations, mostly homogeneous and steady. Most computations involve velocity and position gradient tensors (VGTs and PGTs, respectively). In the context of homogeneous deformations, PGTs are also called finite deformation tensors. This library also offers dynamics for rigid (Jeffery, 1922) and deformable (Eshelby, 1957; Bilby et al., 1975) ellipsoids in slow viscous flows.



### MISCELLANY ###

#' Position gradient tensor for a constant velocity gradient tensor.
#' 
#' Given a VGT L, returns the corresponding PGT F = exp L. This F represents the finite deformation at time t = 1 of a steady, homogeneous progressive deformation that started at t = 0. Requires R package 'expm'.
defExp <- expm

#' Constant velocity gradient tensor deduced from a position gradient tensor.
#' 
#' Given a PGT F, returns the 'principal' logarithm L = log F. This L is a VGT representing a steady, homogeneous progressive deformation, whose finite deformation from t = 0 to t = 1 is F. Requires R package 'expm'.
defLog <- logm

#' Kinematic vorticity of a steady, homogeneous deformation.
#' 
#' The kinematic vorticity, commonly denoted w_k, is the ratio of rotation to distortion in the deformation. If there is no distortion, then it is undefined. Otherwise, it can take on any non-negative real value. If the VGT has only real eigenvalues, then w_k <= 1.
#' @param vgt A real 3x3 matrix. The VGT of the deformation.
#' @return A real number (non-negative) or Inf.
defKinematicVorticity <- function(vgt) {
  d <- 0.5 * (vgt + t(vgt))
  w <- 0.5 * (vgt - t(vgt))
  dNorm <- sqrt(tr(crossprod(d, d)))
  wNorm <- sqrt(tr(crossprod(w, w)))
  wNorm / dNorm
}

#' Finite strain ellipsoid of a homogeneous deformation.
#' 
#' @param pgt A real 3x3 matrix, with positive determinant. The PGT of the deformation.
#' @param doNormalize Logical. Whether to scale the ellipsoid to have the same volume as the unit sphere, or not.
#' @return An ellipsoid (a list with members $tensor, $vector, $a, $logA, $rotation; see ellipsoid.R).
defFiniteStrainEllipsoid <- function(pgt, doNormalize=TRUE) {
  eig <- eigen(pgt %*% t(pgt), symmetric=TRUE)
  # Form a rotation matrix with the finite strain axes as the rows.
  rot <- t(eig$vectors)
  if (det(rot) < 0)
    rot[3,] <- -rot[3,]
  # The Finger tensor has eigenvalues ai^2.
  logA <- 0.5 * log(eig$values)
  ellEllipsoidFromRotationLogA(rot, logA, doNormalize=doNormalize)
}



### INCLINING DEFORMATIONS ###

defInclinedRotation <- function(strike, dip) {
  sinStr <- sin(strike)
  cosStr <- cos(strike)
  sinDip <- sin(dip)
  cosDip <- cos(dip)
  cbind(
    c(sinStr, cosStr, 0),
    c(-cosStr * sinDip, sinStr * sinDip, -cosDip),
    c(-cosStr * cosDip, sinStr * cosDip, sinDip))
}

defInclinedVector <- function(dip, aoc, mag) {
  mag * c(cos(aoc), -sin(dip) * sin(aoc), -cos(dip) * sin(aoc))
}

# Global coordinates y are related to local coordinates x by y = R x + o, so x = R^T (y - o).
# That is, o is the origin of the x-coordinate frame, rendered in y-coordinates.
defInclinedXFromY <- function(y, r, o) {
  as.numeric(t(r) %*% (y - o))
}

defInclinedYFromX <- function(x, r, o) {
  as.numeric(r %*% x + o)
}

defInclinedTensorInXFromInY <- function(tensor, r) {
  t(r) %*% tensor %*% r
}

defInclinedTensorInYFromInX <- function(tensor, r) {
  r %*% tensor %*% t(r)
}

#' Convenience function giving global VGT or PGT for inclined version of homogeneous transpressions.
#' 
#' @param strike A real number. The strike of the shear plane, in radians, but otherwise as usual in geology: measured clockwise from north.
#' @param dip A real number. The dip of the shear plane, in radians. We assume that strike and dip obey the right-hand rule.
#' @param aoc A real number. The angle of oblique convergence, in radians. Kind of like in the monoclinic case.
#' @param mag A real number. The magnitude. Kind of like in the monoclinic case.
#' @param func An R function of the same interface as defTriclinicVGT, defTriclinicPGT, defTransportVGT, etc., giving the tensor in local coordinates.
#' @param hOf1 A real number, positive. The final half-width of the shear zone.
#' @param yOfT A real 3D vector.
#' @param tt A real number. The time t at which y(t) = yOfT.
#' @param origin A real 3D vector. The shear zone origin o in y = R x + o.
#' @return A real 3x3 matrix. The output of func transformed to global coordinates.
defInclinedTensor <- function(strike, dip, aoc, mag, func, hOf1=1, yOfT=c(0, 0, 0), tt=1, origin=c(0, 0, 0)) {
  rot <- defInclinedRotation(strike, dip)
  vec <- defInclinedVector(dip, aoc, mag)
  xOfT <- defInclinedXFromY(yOfT, rot, origin)
  tensor <- func(vec, hOf1=hOf1, xOfT=xOfT, tt=tt)
  defInclinedTensorInYFromInX(tensor, rot)
}



### HETEROGENEOUS TRANSPRESSION (JAEGER, 1962; ROBIN AND CRUDEN, 1994) ###

#' Heterogeneous triclinic transpression (Robin and Cruden, 1994).
#' 
#' @param u A real 3D vector. The boundary plane at x2 = h moves with velocity vector h u. The boundary plane at x2 = -h moves with velocity -h u.
#' @param hOf1 A real number (positive). The half width of the shear zone at time t = 1.
#' @param xOfT A real 3D vector. The position of a particle at time t = 1.
#' @param tt A real number. The time t at which the velocity gradient tensor is desired.
#' @return A real 3x3 matrix (trace zero). The velocity gradient tensor.
defHeterogeneousVGT <- function(u, hOf1, xOfT, tt) {
  hOfT <- hOf1 * exp(u[[2]] * (tt - 1))
  if (xOfT[[2]]^2 >= hOfT^2)
    diag(c(0, 0, 0))
  else {
    l22 <- 1.5 * u[[2]] * (1 - xOfT[[2]]^2 / hOfT^2)
    l32 <- u[[3]] + 3 * u[[2]] * xOfT[[2]] * xOfT[[3]] / hOfT^2
    rbind(
      c(0, u[[1]], 0),
      c(0, l22, 0),
      c(0, l32, -l22))
  }
}



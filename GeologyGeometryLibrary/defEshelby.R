


# Copyright 2016-2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### DEFORMABLE ELLIPSOIDS ###

# Computes the Eshelby J-tensors, assuming a1 = a2 = a3.
defSphereEshelbyJ <- function(a1) {
  jTwo <- (4 * pi) / (15 * a1^2) * rbind(c(3, 1, 1), c(1, 3, 1), c(1, 1, 3))
  jOne <- (4 * pi / 3) * c(1, 1, 1)
  list(jTwo, jOne)
}

# Computes the Eshelby J-tensors, assuming a1 = a2 > a3.
defOblateEshelbyJ <- function(a1, a3) {
  j1 <- 2 * pi * a1^2 * a3 / (a1^2 - a3^2)^(3 / 2) * (acos(a3 / a1) - a3 / a1 * (1 - a3^2 / a1^2)^(1 / 2))
  j2 <- j1
  j3 <- 4 * pi - j1 - j2
  j13 <- (j3 - j1) / (3 * (a1^2 - a3^2))
  j12 <- pi / (3 * a1^2) - j13 / 4
  j11 <- 3 * j12
  j23 <- (j3 - j2) / (3 * (a1^2 - a3^2))
  j22 <- 4 * pi / (3 * a1^2) - j12 - j23
  j33 <- 4 * pi / (3 * a3^2) - j13 - j23
  jTwo <- rbind(c(j11, j12, j13), c(j12, j22, j23), c(j13, j23, j33))
  jOne <- c(j1, j2, j3)
  list(jTwo, jOne)
}

# Computes the Eshelby J-tensors, assuming a1 > a2 = a3.
defProlateEshelbyJ <- function(a1, a2) {
  j3 <- 2 * pi * a1 * a2^2 / (a1^2 - a2^2)^(3 / 2) * (a1 / a2 * (a1^2 / a2^2 - 1)^(1 / 2) - acosh(a1 / a2))
  j2 <- j3
  j1 <- 4 * pi - j2 - j3
  j13 <- (j3 - j1) / (3 * (a1^2 - a2^2))
  j12 <- (j2 - j1) / (3 * (a1^2 - a2^2))
  j11 <- 4 * pi / (3 * a1^2) - j12 - j13
  j23 <- pi / (3 * a2^2) - j12 / 4
  j22 <- 3 * j23
  j33 <- 4 * pi / (3 * a2^2) - j13 - j23
  jTwo <- rbind(c(j11, j12, j13), c(j12, j22, j23), c(j13, j23, j33))
  jOne <- c(j1, j2, j3)
  list(jTwo, jOne)
}

# Helper function for defTypicalEshelbyJOne. Not my best work.
defTypicalEshelbyIntegral <- function(ai, a1, a2, a3) {
  f <- function(u) {
    1 / ((ai^2 + u) * ((a1^2 + u) * (a2^2 + u) * (a3^2 + u))^(1 / 2))
  }
  tryCatch(
    integrate(f, 0, Inf)$value,
    error=function(e) tryCatch(
      integrate(f, 0, 1)$value + integrate(f, 1, Inf)$value,
      error=function(e) tryCatch(
        integrate(f, 0, 1)$value + integrate(f, 1, 10) + integrate(f, 10, Inf)$value,
        error=function(e) {print(paste("error: defTypicalEshelbyIntegral failed on", ai, a1, a2, a3)); NaN})))
}

# Computes the Eshelby J-tensor [J_i], assuming a1 != a2 != a3 != a1. Does not assume a1 > a2 > a3.
defTypicalEshelbyJOne <- function(a) {
  # Set a1 < a2 < a3.
  a1 <- min(a)
  a3 <- max(a)
  a2 <- sum(a) - a1 - a3
  # The integrals for the two largest ai are most likely to be stable?
  j3 <- 2 * pi * a1 * a2 * a3 * defTypicalEshelbyIntegral(a3, a1, a2, a3)
  j2 <- 2 * pi * a1 * a2 * a3 * defTypicalEshelbyIntegral(a2, a1, a2, a3)
  j1 <- 4 * pi - j3 - j2
  # Permute back.
  c(j1, j2, j3)[order(order(a))]
}

# Computes the Eshelby J-tensor [J_ij], assuming a1 != a2 != a3 != a1. Does not assume a1 > a2 > a3.
defTypicalEshelbyJTwo <- function(a, jOne) {
  j13 <- (jOne[[3]] - jOne[[1]]) / (3 * (a[[1]]^2 - a[[3]]^2))
  j12 <- (jOne[[2]] - jOne[[1]]) / (3 * (a[[1]]^2 - a[[2]]^2))
  j11 <- 4 * pi / (3 * a[[1]]^2) - j12 - j13
  j23 <- (jOne[[3]] - jOne[[2]]) / (3 * (a[[2]]^2 - a[[3]]^2))
  j22 <- 4 * pi / (3 * a[[2]]^2) - j12 - j23
  j33 <- 4 * pi / (3 * a[[3]]^2) - j13 - j23
  rbind(c(j11, j12, j13), c(j12, j22, j23), c(j13, j23, j33))
}

# Computes the Eshelby J-tensors, using helper functions for all subcases.
defEshelbyJ <- function(a) {
  if (a[[1]] != a[[2]] && a[[2]] != a[[3]] && a[[3]] != a[[1]]) {
    jOne <- defTypicalEshelbyJOne(a)
    jTwo <- defTypicalEshelbyJTwo(a, jOne)
    list(jTwo, jOne)
  } else if (a[[1]] == a[[2]] && a[[2]] == a[[3]])
    defSphereEshelbyJ(a[[1]])
  else if (a[[1]] == a[[2]] && a[[2]] > a[[3]])
    defOblateEshelbyJ(a[[1]], a[[3]])
  else if (a[[1]] > a[[2]] && a[[2]] == a[[3]])
    defProlateEshelbyJ(a[[1]], a[[3]])
  else if (a[[2]] == a[[3]] && a[[3]] > a[[1]]) {
    js <- defOblateEshelbyJ(a[[2]], a[[1]])
    p <- rbind(c(0, 0, 1), c(0, 1, 0), c(1, 0, 0))
    list(p %*% js[[1]] %*% p, p %*% js[[2]])
  } else if (a[[2]] > a[[3]] && a[[3]] == a[[1]]) {
    js <- defProlateEshelbyJ(a[[2]], a[[1]])
    p <- rbind(c(0, 1, 0), c(1, 0, 0), c(0, 0, 1))
    list(p %*% js[[1]] %*% p, p %*% js[[2]])
  } else if (a[[3]] == a[[1]] && a[[1]] > a[[2]]) {
    js <- defOblateEshelbyJ(a[[3]], a[[2]])
    p <- rbind(c(1, 0, 0), c(0, 0, 1), c(0, 1, 0))
    list(p %*% js[[1]] %*% p, p %*% js[[2]])
  } else { # a[[3]] > a[[1]] && a[[1]] == a[[2]]
    js <- defProlateEshelbyJ(a[[3]], a[[1]])
    p <- rbind(c(0, 0, 1), c(0, 1, 0), c(1, 0, 0))
    list(p %*% js[[1]] %*% p, p %*% js[[2]])
  }
}

#' Deformable ellipsoids of Eshelby (1957) and Bilby et al. (1975).
#' 
#' @param q A real 3x3 matrix (special orthogonal). The orientation of the ellipsoidal inclusion. That is, the rows of Q are unit vectors pointing along the semi-axes of the ellipsoid in a right-handed manner.
#' @param a A real 3D vector (all entries positive). The magnitude of the ellipsoidal inclusion. That is, these are the semi-axis lengths corresponding to the directions in Q.
#' @param d A real 3x3 matrix (symmetric). The stretching tensor of the ambient deformation.
#' @param w A real 3x3 matrix (anti-symmetric). The vorticity tensor of the ambient deformation.
#' @param r A real number (positive). The viscosity ratio between the inclusion and the surrounding material. For example, a passive inclusion has r = 1 and a competent inclusion might have r = 10 or r = 100.
#' @return A real 3x3 matrix. The velocity gradient tensor of the deforming ellipsoidal inclusion.
defEshelbyVGT <- function(q, a, d, w, r) {
  js <- defEshelbyJ(a)
  jTwo <- js[[1]]
  jOne <- js[[2]]
  # Get the diagonal entries of C-tilde into cc.
  dTilde <- q %*% d %*% t(q)
  aa <- diag(c(1, 1, 1)) + (3 / (4 * pi)) * (r - 1) * jTwo %*% diag(a^2)
  bb <- diag(dTilde)
  cc <- solve(aa, bb)
  # Compute the other entries of C-tilde.
  cc12 <- dTilde[[1, 2]] / (1 + (r - 1) * (3 / (4 * pi)) * (a[[1]]^2 + a[[2]]^2) * jTwo[[1, 2]])
  cc13 <- dTilde[[1, 3]] / (1 + (r - 1) * (3 / (4 * pi)) * (a[[1]]^2 + a[[3]]^2) * jTwo[[1, 3]])
  cc23 <- dTilde[[2, 3]] / (1 + (r - 1) * (3 / (4 * pi)) * (a[[2]]^2 + a[[3]]^2) * jTwo[[2, 3]])
  cTilde <- rbind(c(cc[[1]], cc12, cc13), c(cc12, cc[[2]], cc23), c(cc13, cc23, cc[[3]]))
  # Compute K-tilde - W-tilde without computing either of those.
  kTildeMinusWTildeMinusCTilde12 <- (r - 1) / (4 * pi) * (jOne[[1]] - jOne[[2]]) * cTilde[[1, 2]]
  kTildeMinusWTildeMinusCTilde13 <- (r - 1) / (4 * pi) * (jOne[[1]] - jOne[[3]]) * cTilde[[1, 3]]
  kTildeMinusWTildeMinusCTilde23 <- (r - 1) / (4 * pi) * (jOne[[2]] - jOne[[3]]) * cTilde[[2, 3]]
  kTildeMinusWTilde <- cTilde + rbind(
    c(0, kTildeMinusWTildeMinusCTilde12, kTildeMinusWTildeMinusCTilde13),
    c(-kTildeMinusWTildeMinusCTilde12, 0, kTildeMinusWTildeMinusCTilde23),
    c(-kTildeMinusWTildeMinusCTilde13, -kTildeMinusWTildeMinusCTilde23, 0))
  # Return K = Q^T (K-tilde - W-tilde) Q + W.
  t(q) %*% kTildeMinusWTilde %*% q + w
}

#' Simulating the deformation of a deformable ellipsoid in a slow viscous flow.
#' 
#' See Eshelby, (1957); Bilby et al. (1975). Uses the classical fourth-order Runge-Kutta algorithm, based on the differential equation E-dot = -K^T E - E K (Davis et al., 2013).
#' @param e0 A real 3x3 matrix (symmetric, positive-definite). The ellipsoid tensor E at time t = 0.
#' @param l A real 3x3 matrix. Usually trace-zero. The VGT of the ambient deformation.
#' @param r A real number (positive). The viscosity ratio between the clast and the host rock. For example, a competent clast might have r = 10.
#' @param n A real number (positive integer). The number of time steps to use in the simulation.
#' @return A real 3x3 matrix (symmetric, positive-definite). The ellipsoid tensor E at time t = 1.
defClassicalEshelby <- function(e0, l, r, n) {
  d <- 0.5 * (l + t(l))
  w <- 0.5 * (l - t(l))
  eDot <- function(s, e) {
    rotA <- ellRotationAFromTensor(e)
    q <- defDynamicsCorrection(rotA$rotation, rotA$a, d)
    k <- defEshelbyVGT(q, rotA$a, d, w, r)
    -t(k) %*% e - e %*% k
  }
  rungeKutta(e0, eDot, n)
}

#' Simulating the deformation of a deformable ellipsoid in a slow viscous flow.
#' 
#' See Eshelby, (1957); Bilby et al. (1975). Similar to defClassicalEshelby, but uses the left-invariant Lie group, fourth-order Runge-Kutta method of Davis et al. (2013). Should be faster (more precise per time required). Requires R package 'expm'. Validated against my Mathematica code on 2017/07/30.
#' @param e0 A real 3x3 matrix (symmetric, positive-definite). The ellipsoid tensor E at time t = 0.
#' @param l A real 3x3 matrix. Usually trace-zero. The VGT of the ambient deformation.
#' @param r A real number (positive). The viscosity ratio between the clast and the host rock. For example, a competent clast might have r = 10.
#' @param n A real number (positive integer). The number of time steps to use in the simulation.
#' @return A real 3x3 matrix (symmetric, positive-definite). The ellipsoid tensor E at time t = 1.
defLeftEshelby <- function(e0, l, r, n) {
  d <- 0.5 * (l + t(l))
  w <- 0.5 * (l - t(l))
  vel <- function(fInv) {
    rotA <- ellRotationAFromTensor(t(fInv) %*% e0 %*% fInv)
    q <- defDynamicsCorrection(rotA$rotation, rotA$a, d)
    k <- defEshelbyVGT(q, rotA$a, d, w, r)
    -k
  }
  fnInv <- rungeKuttaLeft(diag(c(1, 1, 1)), vel, n)
  t(fnInv) %*% e0 %*% fnInv
}

# Tests.
# ellClassicalEshelby(diag(c(1, 1, 1)), rbind(c(1, 5, -2), c(1, 1, 0), c(0, -3, -1)), 10, 4)
# ellLeftEshelby(diag(c(1, 1, 1)), rbind(c(1, 5, -2), c(1, 1, 0), c(0, -3, -1)), 10, 4)



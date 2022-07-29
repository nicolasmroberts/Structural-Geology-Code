


# Copyright 2016-2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### RIGID ELLIPSOIDS ###

# When the ellipsoid is a spheroid, its orientation is ambiguous. However, the dynamical calculations assume that the orientation has been chosen to line up with the eigenvectors of D in a particular way. So this function returns the adjusted orientation. In practice this adjustment is important mainly when using rigid spheroids or when starting a deformable simulation from a spheroid.
defDynamicsCorrection <- function(q, a, d) {
  qNew <- q
  if (a[[1]] == a[[2]]) {
    b <- qNew %*% d %*% t(qNew)
    if (b[[2, 2]] == b[[1, 1]])
      s <- pi / 4
    else
      s <- 0.5 * atan(2 * b[[1, 2]] / (b[[2, 2]] - b[[1, 1]]))
    r <- matrix(c(cos(s), sin(s), 0, -sin(s), cos(s), 0, 0, 0, 1), 3, 3)
    qNew <- r %*% qNew
  }
  if (a[[1]] == a[[3]]) {
    b <- qNew %*% d %*% t(qNew)
    if (b[[1, 1]] == b[[3, 3]])
      s <- pi / 4
    else
      s <- 0.5 * atan(2 * b[[1, 3]] / (b[[1, 1]] - b[[3, 3]]))
    r <- matrix(c(cos(s), 0, -sin(s), 0, 1, 0, sin(s), 0, cos(s)), 3, 3)
    qNew <- r %*% qNew
  }
  if (a[[2]] == a[[3]]) {
    b <- qNew %*% d %*% t(qNew)
    if (b[[3, 3]] == b[[2, 2]])
      s <- pi / 4
    else
      s <- 0.5 * atan(2 * b[[2, 3]] / (b[[3, 3]] - b[[2, 2]]))
    r <- matrix(c(1, 0, 0, 0, cos(s), sin(s), 0, -sin(s), cos(s)), 3, 3)
    qNew <- r %*% qNew
  }
  qNew
}

# The velocity gradient tensor K = V of the Jeffery (1922) rigid ellipsoid. This is not the same thing as Qdot.
defJefferyVGT <- function(q, a, d, w) {
  q <- defDynamicsCorrection(q, a, d)
  dTilde <- q %*% d %*% t(q)
  wTilde12 <- (a[[1]]^2 - a[[2]]^2) / (a[[1]]^2 + a[[2]]^2) * dTilde[[1, 2]]
  wTilde13 <- (a[[1]]^2 - a[[3]]^2) / (a[[1]]^2 + a[[3]]^2) * dTilde[[1, 3]]
  wTilde23 <- (a[[2]]^2 - a[[3]]^2) / (a[[2]]^2 + a[[3]]^2) * dTilde[[2, 3]]
  wTilde <- matrix(c(0, -wTilde12, -wTilde13, wTilde12, 0, -wTilde23, wTilde13, wTilde23, 0), 3, 3)
  w - t(q) %*% wTilde %*% q
}

#' Simulating the rotation of a rigid ellipsoid in a slow viscous flow.
#' 
#' See Jeffery (1922). Uses the classical fourth-order Runge-Kutta algorithm.
#' @param q0 A real 3x3 matrix (special orthogonal). The ellipsoid's orientation Q at time t = 0.
#' @param a A real 3-dimensional vector. The semi-axis lengths a1, a2, a3, in order corresponding to the rows of Q.
#' @param l A real 3x3 matrix. Usually trace-zero. The VGT of the ambient deformation.
#' @param n A real number (positive integer). The number of time steps to use in the simulation.
#' @return A real 3x3 matrix (special orthogonal). The ellipsoid's orientation Q at time t = 1.
defClassicalJeffery <- function(q0, a, l, n) {
  d <- 0.5 * (l + t(l))
  w <- 0.5 * (l - t(l))
  vel <- function(s, q) {-q %*% defJefferyVGT(q, a, d, w)}
  rungeKutta(q0, vel, n, rotProjectedMatrix)
}

#' Simulating the rotation of a rigid ellipsoid in a slow viscous flow (Jeffery, 1922).
#' 
#' The same calculation as defClassicalJeffery, but using a left-invariant Lie group Runge-Kutta method instead of a classical Runge-Kutta method. An improved version of what's in Davis et al. (2013). This function may deliver less error (per time spent) than defClassicalJeffery, but I have not done testing to confirm that.
#' @param q0 A real 3x3 matrix (special orthogonal). The ellipsoid's orientation Q at time t = 0.
#' @param a A real 3-dimensional vector. The semi-axis lengths a1, a2, a3, in order corresponding to the rows of Q.
#' @param l A real 3x3 matrix. Usually trace-zero. The VGT of the ambient deformation.
#' @param n A real number (positive integer). The number of time steps to use in the simulation.
#' @return A real 3x3 matrix (special orthogonal). The ellipsoid's orientation Q at time t = 1.
defLeftJeffery <- function(q0, a, l, n) {
  d <- 0.5 * (l + t(l))
  w <- 0.5 * (l - t(l))
  vel <- function(q) {-defJefferyVGT(q, a, d, w)}
  rungeKuttaLeft(q0, vel, n, exponential=rotExp)
}

#' Simulating the rotation of a rigid ellipsoid in a slow viscous flow (Jeffery, 1922).
#' 
#' The same calculation as defClassicalJeffery, but using a right-invariant Lie group Runge-Kutta method instead of a classical Runge-Kutta method. An improved version of what's in Davis et al. (2013). This function may deliver less error (per time spent) than defClassicalJeffery, but I have not done testing to confirm that.
#' @param q0 A real 3x3 matrix (special orthogonal). The ellipsoid's orientation Q at time t = 0.
#' @param a A real 3-dimensional vector. The semi-axis lengths a1, a2, a3, in order corresponding to the rows of Q.
#' @param l A real 3x3 matrix. Usually trace-zero. The VGT of the ambient deformation.
#' @param n A real number (positive integer). The number of time steps to use in the simulation.
#' @return A real 3x3 matrix (special orthogonal). The ellipsoid's orientation Q at time t = 1.
defRightJeffery <- function(q0, a, l, n) {
  d <- 0.5 * (l + t(l))
  w <- 0.5 * (l - t(l))
  vel <- function(s, qT) {defJefferyVGT(t(qT), a, d, w)}
  t(rungeKuttaRight(t(q0), vel, n, exponential=rotExp))
}

# Test that all three functions work identically or close to it.
#q0 <- rotUniform()
#a <- abs(rnorm(3))
#a <- a / prod(a)^(1 / 3)
#l <- replicate(3, rnorm(3))
#l <- l - diag(c(1, 1, 1)) * tr(l) / 3
#defClassicalJeffery(q0, a, l, 40)
#defLeftJeffery(q0, a, l, 10)
#defRightJeffery(q0, a, l, 10)



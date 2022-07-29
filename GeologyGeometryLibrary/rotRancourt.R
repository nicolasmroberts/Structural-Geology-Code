


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### RANCOURT ET AL. (2000) ###

#' Approximation of one rotation in the tangent space at the other rotation.
#'
#' This function is based on Rancourt et al. (2000). It produces a somewhat different approximation from that of rotLeftTangentFromMatrix. Appropriate only if r and center are close to each other. Definitely inappropriate if the distance between r and center is greater than pi / 2.
#' @param r A rotation matrix.
#' @param center A rotation matrix.
#' @return A 3D real vector.
rotRancourtFromMatrix <- function(r, center) {
  w <- (crossprod(center, r) - crossprod(r, center)) / 2.0
  rotVectorFromAntisymmetric(w)
}

#' Mapping a tangent space into the space of rotations.
#'
#' Inverse to rotRancourtFromMatrix.
#' @param v A 3D real vector.
#' @param center A rotation matrix.
#' @return A rotation matrix.
rotMatrixFromRancourt <- function(v, center) {
  cosine <- squareRoot(1 - dot(v, v))
  w <- rotAntisymmetricFromVector(v)
  r <- diag(c(1, 1, 1)) + w + (w %*% w) / (1 + cosine)
  center %*% r
}

#' P-value function for any hypothesized mean, using tangent space approximation of Rancourt et al. (2000).
#' 
#' Appropriate only if the rotations are tightly clustered. Definitely inappropriate if any rotation is more than pi / 2 away from the mean. Uses Eq. (2.6) of Rancourt et al. (2000), 'Using orientation statistics to investigate variations in human kinematics'. Inserts an extra factor of n on the left side of that equation. Without this extra factor, coverage rates are much too large. Also the corresponding author confirmed via e-mail on 2015/08/31 that the extra factor should be there.
#' @param rs A list of rotation matrices. The data.
#' @return An R function from {rotation matrices} to {real numbers union NA}. For any given hypothesized mean, returns NA if that mean is more than pi / 2 away from the mean of the rs, or the p-value if not.
rotRancourtInference <- function(rs) {
  udv <- rotSVD(arithmeticMean(rs))
  mHat <- udv$u %*% t(udv$v)
  vs <- lapply(rs, rotRancourtFromMatrix, mHat)
  ms <- lapply(vs, function(v) {outer(v, v)})
  n <- length(rs)
  s <- arithmeticMean(ms) * n / (n - 1)
  sInv <- solve(s)
  # This line has an extra n factor, compared to Eq. (2.6) of Rancourt et al. (2000).
  g <- n * (n - 3) / (3 * (n - 1)) * sInv
  f <- function(r) {
    if (rotDistance(r, mHat) > pi / 2)
      NA
    else {
      v <- rotRancourtFromMatrix(r, mHat)
      1 - pf(v %*% g %*% v, 3, n - 3)
    }
  }
  f
}



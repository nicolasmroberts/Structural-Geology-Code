


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### UNIFORM DISTRIBUTION ###

#' Sampling from the uniform distribution on the space of rotations.
#'
#' @param n A real number (positive integer) or NULL.
#' @return Either a rotation matrix (if n is NULL) or a list of n rotation matrices (if n is a number).
rotUniform <- function(n=NULL) {
  if (is.null(n)) {
    x <- rayNormalized(c(1, 0, 0) - rayUniform())
    refl <- diag(c(1, 1, 1)) - 2 * (x %o% x)
    a <- runif(1, 0, 2 * pi)
    -refl %*% rotMatrixAboutX(a)
  }
  else
    replicate(n, rotUniform(), simplify=FALSE)
}

#' Test of uniformity, based on the Rayleigh statistic.
#'
#' @param rs A list of rotation matrices.
#' @return A list with members $p, $rayleigh. $p is a real number (between 0 and 1), the p-value for the test. Low values indicate that the sample was not drawn from the uniform distribution. $rayleigh is the Rayleigh statistic that produced $p.
rotRayleighInference <- function(rs) {
  rBar <- arithmeticMean(rs)
  rayleigh <- 3 * length(rs) * tr(crossprod(rBar, rBar))
  p <- 1 - pchisq(rayleigh, 9)
  list(p=p, rayleigh=rayleigh)
}



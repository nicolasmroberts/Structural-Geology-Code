


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# A line is expressed as a unit 3D vector in Cartesian coordinates, as an R vector u = c(x, y, z) where x^2 + y^2 + z^2 == 1. It is crucial to remember that u and -u represent the same line. Many of the line functions here are implemented in terms of the ray functions of rays.R. In particular, when lines are tightly concentrated, you can often ignore their 'negative copies' and treat them as rays, with no appreciable effect on the statistics.



### UNIFORM DISTRIBUTION ###

#' Uniformly random lines.
#' 
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single line. If n is a positive integer, then a list of n lines.
lineUniform <- function(n=NULL) {
  if (is.null(n))
    lower(rayUniform())
  else
    lapply(rayUniform(n), lower)
}

#' Test of uniformity based on the modified Bingham statistic.
#' 
#' See Mardia and Jupp (2000, Section 10.7.1). Not consistent against alternatives with scatter matrix Tbar = I / dim, where dim is the ambient dimension (i.e. data are points in RP^(dim - 1)).
#' @param ls A list of lines.
#' @return A real number. The p-value for the null hypothesis that the lines are uniformly distributed.
lineUniformBinghamPValue <- function(ls) {
  n <- length(ls)
  p <- length(ls[[1]])
  # Compute the Bingham statistic S.
  scat <- lineMeanScatter(ls)$values
  trTBarSq <- sum(scat^2)
  s <- (trTBarSq - 1 / p) * n * p * (p + 2) / 2
  # Compute the modified Bingham statistic S*.
  b0 <- (2 * p^2 + 3 * p + 4) / (6 * (p + 4))
  b1 <- -(4 * p^2 + 3 * p - 4) / (3 * (p + 4) * (p^2 + p + 2))
  b2 <- -4 * (p^2 - 4) / (3 * (p + 4) * (p^2 + p + 2) * (p^2 + p + 6))
  sStar <- s * (1 - (b0 + b1 * s + b2 * s^2) / n)
  # Then S* is distributed like chi squared with (p - 1) (p + 1) / 2 degrees of freedom.
  pchisq(sStar, (p - 1) * (p + 1) / 2, lower.tail=FALSE)
}

#lineUniformBinghamPValue(lineUniform(100))
#lineUniformBinghamPValue(lineWatson(mu=lineUniform(), kappa=4, n=100))

# To do: Gine test of uniformity, Mardia and Jupp (2000, p. 233)



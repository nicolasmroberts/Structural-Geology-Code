


# Copyright 2016-2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### COMPUTING SPO ###

# poleRakeOther is a rotation matrix with rows pointing along pole, rake, and other direction.
# Returns coefficients of B11, B12, B13, B23, B22, 1, and extra variable C in three equations.
ellRobinCoefficients <- function(poleRakeOther, rakeSemiaxisLength, otherSemiaxisLength) {
  # l is the matrix with columns pole, rake, other. Only its second and third columns will be used.
  l <- t(poleRakeOther)
  # First equation, based on (1, 1) entry.
  first <- c(
    l[[1, 2]]^2 - l[[3, 2]]^2,
    2 * l[[1, 2]] * l[[2, 2]],
    2 * l[[1, 2]] * l[[3, 2]],
    2 * l[[2, 2]] * l[[3, 2]],
    l[[2, 2]]^2 - l[[3, 2]]^2,
    0,
    -rakeSemiaxisLength^-2)
  # Second equation, based on (2, 2) entry.
  second <- c(
    l[[1, 3]]^2 - l[[3, 3]]^2,
    2 * l[[1, 3]] * l[[2, 3]],
    2 * l[[1, 3]] * l[[3, 3]],
    2 * l[[2, 3]] * l[[3, 3]],
    l[[2, 3]]^2 - l[[3, 3]]^2,
    0,
    -otherSemiaxisLength^-2)
  # Third equation, based on (1, 2) or (2, 1) entry.
  third <- c(
    l[[1, 2]] * l[[1, 3]] - l[[3, 2]] * l[[3, 3]],
    l[[1, 3]] * l[[2, 2]] + l[[1, 2]] * l[[2, 3]],
    l[[1, 3]] * l[[3, 2]] + l[[1, 2]] * l[[3, 3]],
    l[[2, 3]] * l[[3, 2]] + l[[2, 2]] * l[[3, 3]],
    l[[2, 2]] * l[[2, 3]] - l[[3, 2]] * l[[3, 3]],
    3 * l[[3, 2]] * l[[3, 3]],
    0)
  list(first, second, third)
}

#' Fit an ellipsoid to elliptical sections, using their shape but not size.
#' 
#' Warning: This function is not well tested. Actually I have reason to believe that it is quite wrong. Anyway, this is the second case treated by Robin (2002). The output ellipsoid tensor is normalized to have trace 3, and might not actually be positive-definite at all.
#' @param poleRakeOthers A list of real 3x3 matrices (special orthogonal). Each matrix describes the orientation of an ellipse in space. The first row is the pole to the plane. The second row is the rake of one of the ellipse's axes in that plane. The third row is the cross product of the first two.
#' @param rakeSemiaxisLengths A vector of real numbers. The length of the semi-axis indicated by the rake in the first argument.
#' @param otherSemiaxisLengths A vector of real numbers. The length of the semi-axis perpendicular to the rake.
#' @return A real 3x3 matrix (symmetric, trace-3). The putative ellipsoid tensor.
ellRobin <- function(poleRakeOthers, rakeSemiaxisLengths, otherSemiaxisLengths) {
  # Construct a system X B = Y of linear equations.
  n <- length(poleRakeOthers)
  x <- matrix(0, 3 * n, 5 + n)
  y <- replicate(3 * n, 0)
  for (i in 1:n) {
    eqns <- ellRobinCoefficients(poleRakeOthers[[i]], rakeSemiaxisLengths[[i]], otherSemiaxisLengths[[i]])
    x[(i * 3 - 2),1:5] <- eqns[[1]][1:5]
    x[(i * 3 - 2),(5 + i)] <- eqns[[1]][7]
    x[(i * 3 - 1),1:5] <- eqns[[2]][1:5]
    x[(i * 3 - 1),(5 + i)] <- eqns[[2]][7]
    x[(i * 3),1:5] <- eqns[[3]][1:5]
    x[(i * 3),(5 + i)] <- eqns[[3]][7]
    y[[i * 3 - 2]] <- -eqns[[1]][[6]]
    y[[i * 3 - 1]] <- -eqns[[2]][[6]]
    y[[i * 3]] <- -eqns[[3]][[6]]
  }
  # Solve for B.
  fit <- lm.fit(x, y)
  # For now, just rebuild the trace-3 ellipsoid tensor.
  es <- fit$coefficients
  rbind(c(es[[1]], es[[2]], es[[3]]),
        c(es[[2]], es[[4]], es[[5]]),
        c(es[[3]], es[[5]], 3 - es[[1]] - es[[4]]))
}

#' Fit an ellipsoid to elliptical sections, using their shape but not size.
#' 
#' This is my custom method for fitting SPO. Unlike the method of Robin (2002, Case 2), this method is guaranteed to produce a positive-definite ellipsoid tensor E. Currently I force volume normalization (det E = 1) as well. Works well, except when minEigenvalue is negative. Works less well if BFGS is replaced with the default (Nelder-Mead).
#' @param poleRakeOthers A list of real 3x3 matrices (special orthogonal). Each matrix describes the orientation of an ellipse in space. The first row is the pole to the plane. The second row is the rake of one of the ellipse's axes in that plane. The third row is the cross product of the first two.
#' @param rakeSemiaxisLengths A vector of real numbers. The length of the semi-axis indicated by the rake in the first argument.
#' @param otherSemiaxisLengths A vector of real numbers. The length of the semi-axis perpendicular to the rake.
#' @param numSteps A real number (positive integer). The number of steps to use in the optimization algorithm.
#' @return A list with members $ellipsoid, $error, $minEigenvalue, and $value. $ellipsoid is an ellipsoid. $error is an error code; if it is non-zero, then an error occurred; try increasing numSteps. $minEigenvalue is the least eigenvalue of the Hessian at the putative optimum; if it is non-positive, then an error occurred. $value is the value of the misfit function at the optimum.
ellSPO <- function(poleRakeOthers, rakeSemiaxisLengths, otherSemiaxisLengths, numSteps=10000) {
  # Pre-process the data.
  n <- length(poleRakeOthers)
  ls <- lapply(poleRakeOthers, function(r) r[2:3,])
  bs <- thread(function(a1, a2) diag(c(a1, a2)^-2), rakeSemiaxisLengths, otherSemiaxisLengths)
  # Build the misfit function to be minimized.
  misfit <- function(pars) {
    s <- rbind(
      c(pars[[1]], pars[[2]], pars[[3]]), 
      c(pars[[2]], pars[[4]], pars[[5]]), 
      c(pars[[3]], pars[[5]], -pars[[1]] - pars[[4]]))
    e <- ellExp(s)
    diffs <- lapply(1:n, function(i) {exp(pars[[5 + i]]) * bs[[i]] - ls[[i]] %*% e %*% t(ls[[i]])})
    normSqs <- sapply(diffs, function(diff) tr(t(diff) %*% diff))
    sum(normSqs)
  }
  # Seed the minimization from the unit sphere and all k_i = 0.5 arbitrarily.
  seed <- c(0, 0, 0, 0, 0, replicate(n, 0.5))
  solution <- optim(seed, misfit, hessian=TRUE, method="BFGS", control=list(maxit=numSteps))
  # Report the answer and diagnostic information.
  eigvals <- eigen(solution$hessian, symmetric=TRUE, only.values=TRUE)$values
  s <- rbind(
    c(solution$par[[1]], solution$par[[2]], solution$par[[3]]), 
    c(solution$par[[2]], solution$par[[4]], solution$par[[5]]), 
    c(solution$par[[3]], solution$par[[5]], -solution$par[[1]] - solution$par[[4]]))
  ell <- ellEllipsoidFromTensor(ellExp(s), doNormalize=TRUE)
  list(ellipsoid=ell, error=solution$convergence, minEigenvalue=min(eigvals), value=solution$value)
}

# Testing for ellSPO. Noiseless. Test exactly mirrors the optimization procedure.
ellSPOTest <- function(n) {
  # Make a random ellipsoid.
  q <- rotUniform()
  a <- exp(rnorm(3))
  e <- t(q) %*% diag(a^-2) %*% q
  # Make n random sections.
  f <- function() {
    l <- rotUniform()[2:3,]
    b <- l %*% e %*% t(l)
    eig <- eigen(b, symmetric=TRUE)
    l <- t(eig$vectors) %*% l
    b <- l %*% e %*% t(l)
    rake <- b[[1, 1]]^(-1 / 2)
    other <- b[[2, 2]]^(-1 / 2)
    list(l=l, rake=rake, other=other)
  }
  sections <- replicate(n, f(), simplify=FALSE)
  # Dissect the sections into the format desired by ellRobin and ellSPO.
  poleRakeOthers <- lapply(sections, 
                           function(s) rbind(cross(s$l[1,], s$l[2,]), s$l[1,], s$l[2,]))
  rakeSemiaxisLengths <- sapply(sections, function(s) s$rake)
  otherSemiaxisLengths <- sapply(sections, function(s) s$other)
  # Compare the true answer to the deduced answer.
  print("true ellipsoid tensor, volume-normalized:")
  print(e * det(e)^(-1 / 3))
  print("ellSPO result:")
  pred <- ellSPO(poleRakeOthers, rakeSemiaxisLengths, otherSemiaxisLengths)
  print(pred$ellipsoid$tensor)
  print(c(pred$error, pred$minEigenvalue))
  print("ellRobin result, volume-normalized:")
  pred <- ellRobin(poleRakeOthers, rakeSemiaxisLengths, otherSemiaxisLengths)
  print(pred * det(pred)^(-1 / 3))
}

# This is an example hand-ported from Mathematica. Again noiseless, but the sections are not being generated by the same code that does the optimization. We should get rbind(c(0.879642, -0.0768036, -0.0419311), c(-0.0768036, 1.06686, 0.0109123), c(-0.0419311, 0.0109123, 1.07437)).
ellSPOTestSpecific <- function() {
  rakeOthers <- list(
    rbind(
      c(-0.336673, 0.941231, 0.0271225),
      c(-0.941158, -0.335463, -0.0410727)),
    rbind(
      c(-0.251698, 0.938829, 0.23505),
      c(-0.506325, 0.0792426, -0.858694)),
    rbind(
      c(-0.263783, 0.947118, -0.182718),
      c(0.865012, 0.31609, 0.389668)),
    rbind(
      c(0.065659, -0.526426, 0.847682),
      c(-0.324202, -0.814681, -0.48082)))
  rakeSemiaxisLengths <- c(0.955355, 0.952817, 0.960189, 0.970211)
  otherSemiaxisLengths <- c(1.08491, 1.00371, 1.07812, 0.998093)
  poleRakeOthers <- lapply(
    rakeOthers, function(ro) rotProjectedMatrix(rbind(cross(ro[1,], ro[2,]), ro[1,], ro[2,])))
  ellSPO(poleRakeOthers, rakeSemiaxisLengths, otherSemiaxisLengths)
}



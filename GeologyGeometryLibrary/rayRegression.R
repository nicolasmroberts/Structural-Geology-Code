


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



#' Fitting a great circle arc to a set of rays, based on an independent scalar variable.
#' 
#' @param xs A vector of real numbers. Values of the independent scalar variable x.
#' @param ls A list of rays, of the same length as xs.
#' @param numSteps A real number (positive integer). A bound on the number of numerical optimization iterations allowed. If you find that the results have $error != 0, then try increasing numSteps.
#' @param numPoints A real number (positive integer). The number of points along the regression curve requested.
#' @return A list (a, rotation, error, minEigenvalue, rSquared, prediction). The great circle curve is v(x) = R [cos(a x) sin(a x) 0]^T, where R == rotation. The third column of R is the pole of the great circle. a is the angle of rotation about that pole per unit x, in radians. So a / degree is the angle in degrees. error should be 0 and minEigenvalue should be positive. Otherwise there was some problem in the optimization. rSquared is the squared correlation coefficient, 0 <= R^2 <= 1. If numPoints > 0 then the return list also has a member $points consisting of numPoints+1 rays. prediction is an R function that takes an x value (real number) as input and produces the ray v(x) as output.
rayGeodesicRegression <- function(xs, ls, numSteps=1000, numPoints=0) {
  # Let l0 be the l whose x is closest to zero.
  l0 <- ls[[which.min(as.numeric(xs)^2)]]
  # Define the function to be minimized.
  n <- length(ls)
  e <- function(wb) {
    a <- rotExp(rotAntisymmetricFromVector(wb[1:3]))
    f <- function(i) {
      rayDistance(ls[[i]], as.numeric(a %*% c(cos(wb[[4]] * xs[[i]]), sin(wb[[4]] * xs[[i]]), 0)))^2
    }
    sum(sapply(1:n, f)) / (2 * n)
  }
  # Find the minimum, using the constant geodesic l0 as the seed.
  seed <- c(0, 0, 0, 0)
  solution <- optim(seed, e, hessian=TRUE, control=list(maxit=numSteps))
  # Report diagnostic information.
  eigvals <- eigen(solution$hessian, symmetric=TRUE, only.values=TRUE)$values
  rot <- rotExp(rotAntisymmetricFromVector(solution$par[1:3]))
  a <- solution$par[[4]]
  rSq <- 1 - solution$value / rayVariance(ls, rayProjectedMean(ls))
  result <- list(a=a, rotation=rot, error=solution$convergence, minEigenvalue=min(eigvals), rSquared=rSq)
  result$prediction <- function(x) {
    as.numeric(rot %*% c(cos(a * x), sin(a * x), 0))
  }
  if (numPoints >= 1)
    result$points <- lapply(seq(from=min(xs), to=max(xs), length.out=(numPoints + 1)), result$prediction)
  result
}

#' Fitting a great circle arc to a set of lines, based on an independent scalar variable.
#' 
#' Similar to rayGeodesicRegression, but internally rescales x to [0, 1] for better performance. Try this function first, and use the other function only if you find reason to.
#' @param xs See rayGeodesicRegression.
#' @param ls See rayGeodesicRegression.
#' @param numSteps See rayGeodesicRegression.
#' @param numPoints See rayGeodesicRegression.
#' @return See rayGeodesicRegression.
rayRescaledGeodesicRegression <- function(xs, ls, numSteps=1000, numPoints=0) {
  x0 <- min(xs)
  x1 <- max(xs)
  # Perform regression in scaled coordinates.
  regr <- rayGeodesicRegression(scales(xs), ls, numSteps, numPoints=0)
  # Scale the results back into the original coordinates.
  a <- regr$a / (x1 - x0)
  col1 <- regr$rotation[,1] * cos(a * x0) - regr$rotation[,2] * sin(a * x0)
  col2 <- regr$rotation[,1] * sin(a * x0) + regr$rotation[,2] * cos(a * x0)
  rot <- rotProjectedMatrix(cbind(col1, col2, cross(col1, col2)))
  result <- list(a=a, rotation=rot, error=regr$error, minEigenvalue=regr$minEigenvalue, rSquared=regr$rSquared)
  result$prediction <- function(x) {
    as.numeric(rot %*% c(cos(a * x), sin(a * x), 0))
  }
  if (numPoints >= 1)
    result$points <- lapply(seq(from=x0, to=x1, length.out=numPoints), result$prediction)
  result
}

# Obsolete; using permutedRSquareds instead.
rayGeodesicRegressionPermutations <- function(xs, ls, numPerms, numSteps=10000) {
  ys <- scales(xs)
  f <- function(i) {
    print(i / numPerms)
    regr <- rayGeodesicRegression(sample(ys, size=length(ys)), ls, numSteps)
    c(regr$error, regr$minEigenvalue, regr$rSquared)
  }
  perms <- sapply(1:numPerms, f)
  perms[3,][perms[1,] == 0 & perms[2,] > 0]
}

# Test.
#r <- rayUniform()
#rs <- lapply(rayGeodesicPoints(r, rayOrthogonalUniform(r), numSteps=10), rayFisher, 20)
#xs <- seq(from=0, to=1, length.out=length(rs))
#regr <- rayRescaledGeodesicRegression(xs, rs, numSteps=1000, numPoints=10)
#regr$error
#regr$minEigenvalue
#regr$rSquared
#regr$a / degree
#geoTrendPlungeDegFromCartesian(regr$rotation[,3])
#rayEqualAreaPlot(rs, colors=hues(xs), curves=list(regr$points))
#perms <- rayGeodesicRegressionPermutations(xs, rs, 1000)
#length(perms)
#sum(perms > regr$rSquared) / length(perms)

# Warning: not well tested. angleBound > 0 is a bound on the amount of rotation (either positive or negative) per unit of x.
raySmallCircleRegression <- function(xs, us, numSeeds=5, numSteps=1000, numPoints=0, angleBound=Inf) {
  f <- function(phiThetaAlphaTauSigma) {
    pole <- cartesianFromSpherical(c(1, phiThetaAlphaTauSigma[1:2]))
    uOf0 <- cartesianFromSpherical(c(1, phiThetaAlphaTauSigma[4:5]))
    pred <- function(x) {
      as.numeric(rotMatrixFromAxisAngle(c(pole, x * phiThetaAlphaTauSigma[[3]])) %*% uOf0)
    }
    preds <- lapply(xs, pred)
    dists <- mapply(rayDistance, preds, us)
    dot(dists, dists) / (2 * length(us))
  }
  # Find the minimum, starting from a few seeds.
  best <- list(value=(2 * pi^2))
  for (i in 1:numSeeds) {
    seed <- c(sphericalFromCartesian(rayUniform())[2:3], runif(1, -pi, pi), sphericalFromCartesian(rayUniform())[2:3])
    if (angleBound == Inf)
      solution <- optim(seed, f, hessian=TRUE, 
                        control=list(maxit=numSteps), method="L-BFGS-B")
    else
      solution <- optim(seed, f, hessian=TRUE, 
                        control=list(maxit=numSteps), method="L-BFGS-B",
                        lower=c(-Inf, -Inf, -angleBound, -Inf, -Inf), 
                        upper=c(Inf, Inf, angleBound, Inf, Inf))
    if (solution$value <= best$value)
      best <- solution
  }
  # Report results.
  eigvals <- eigen(best$hessian, symmetric=TRUE, only.values=TRUE)$values
  pole <- cartesianFromSpherical(c(1, best$par[1:2]))
  angle <- best$par[[3]]
  uOf0 <- cartesianFromSpherical(c(1, best$par[4:5]))
  var <- rayVariance(us, rayProjectedMean(us))
  rSquared <- 1 - best$value / var
  pred <- function(x) {as.numeric(rotMatrixFromAxisAngle(c(pole, x * angle)) %*% uOf0)}
  results <- list(error=best$convergence, minEigenvalue=min(eigvals),
                  pole=pole, angle=angle, rSquared=rSquared, prediction=pred)
  if (numPoints >= 1) {
    ys <- seq(from=min(xs), to=max(xs), length.out=numPoints)
    results$points <- lapply(ys, pred)
  }
  results
}

# Warning: not well tested. angleBound > 0 is a bound on the amount of rotation (either positive or negative) per unit of x.
rayRescaledSmallCircleRegression <- function(xs, ls, angleBound=Inf, ...) {
  # Perform regression in scaled coordinates.
  x0 <- min(xs)
  x1 <- max(xs)
  results <- raySmallCircleRegression(scales(xs), ls, angleBound=(angleBound * (x1 - x0)), ...)
  # Scale the results back into the original coordinates.
  results$rescaledPrediction <- results$prediction
  results$prediction <- function(x) {results$rescaledPrediction((x - x0) / (x1 - x0))}
  results$angle <- results$angle / (x1 - x0)
  results
}

# Obsolete. See permutedRSquareds.
raySmallCircleRegressionPermutations <- function(xs, ls, numPerms, ...) {
  ys <- scales(xs)
  f <- function(i) {
    print(i / numPerms)
    regr <- raySmallCircleRegression(sample(ys, size=length(ys)), ls, ...)
    c(regr$error, regr$minEigenvalue, regr$rSquared)
  }
  perms <- sapply(1:numPerms, f)
  perms[3,][perms[1,] == 0 & perms[2,] > 0]
}

# Demonstration and testing. Permutations don't work very well yet, for some reason.
raySmallCircleRegressionTest <- function(n=10, kappa=100, numPerms=100, numSteps=10000) {
  pole <- rayUniform()
  angle <- runif(1, 0, pi)
  uOf0 <- rayUniform()
  xs <- runif(n, -1, 1)
  xs <- sort(xs)
  pred <- function(x) as.numeric(rotMatrixFromAxisAngle(c(pole, x * angle)) %*% uOf0)
  preds <- lapply(xs, pred)
  #rayEqualAreaPlot(preds, colors=hues(xs))
  us <- lapply(preds, rayFisher, kappa)
  #rayEqualAreaPlot(us, colors=hues(xs))
  #print(xs)
  regr <- raySmallCircleRegression(xs, us, numPoints=20)
  print("error, minEigenvalue, R^2, angle:")
  print(c(regr$error, regr$minEigenvalue, regr$rSquared, regr$angle))
  print("pole:")
  print(regr$pole)
  print(geoTrendPlungeDegFromCartesian(regr$pole))
  rayEqualAreaPlot(us, colors=hues(xs), curves=list(regr$points))
  regr <- rayRescaledSmallCircleRegression(xs, us, numPoints=20)
  print("error, minEigenvalue, R^2, angle:")
  print(c(regr$error, regr$minEigenvalue, regr$rSquared, regr$angle))
  print("pole:")
  print(regr$pole)
  print(geoTrendPlungeDegFromCartesian(regr$pole))
  rayEqualAreaPlot(us, colors=hues(xs), curves=list(regr$points))
  rSquareds <- raySmallCircleRegressionPermutations(xs, us, numPerms=numPerms, numSteps=numSteps)
  print("numPerms, num successfully completed, num with greater R^2:")
  print(c(numPerms, length(rSquareds), sum(rSquareds > regr$rSquared)))
}
#raySmallCircleRegressionTest()






# Copyright 2016-2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# An ellipsoid is a complete description of an ellipsoid in three-dimensional space: its size, shape, and orientation (but not its location). There are three degrees of freedom in orientation and three degrees of freedom in size-shape. Sometimes we normalize ellipsoids to have a particular volume, in which case only two degrees of freedom remain in shape.

# An ellipsoid can be described as a symmetric, positive-definite ellipsoid tensor E, in that the (boundary of the) ellipsoid is the set of points x such that x^T E x == 1. Or it can be described as an orientation and semi-axis lengths, or as a log-ellipsoid vector. In this R code, an 'ellipsoid' is a list of five elements: $vector, $tensor, $a, $logA, $rotation. $a are the semi-axis lengths, and $logA is their logarithms. $rotation is a special orthogonal matrix with the semi-axes of the ellipsoid along its rows. See also geoDataFromFile.

# Because they inhabit a vector space, we can throw all of multivariate statistics at the ellipsoid vectors. For ideas, see http://cran.r-project.org/web/views/Multivariate.html.



### CONVERSIONS AMONG REPRESENTATIONS ###

#' A log-ellipsoid tensor from an ellipsoid tensor, respecting normalization.
#' 
#' @param ell A real 3x3 matrix (symmetric, positive-definite).
#' @param A real 3x3 matrix (symmetric).
ellLog <- function(ell) {
  eig <- eigen(ell, symmetric=TRUE)
  eig$vectors %*% diag(log(eig$values)) %*% t(eig$vectors)
}

#' An ellipsoid tensor from a log-ellipsoid tensor, respecting normalization.
#' 
#' @param logEll A real 3x3 matrix (symmetric).
#' @param A real 3x3 matrix (symmetric, positive-definite).
ellExp <- function(logEll) {
  eig <- eigen(logEll, symmetric=TRUE)
  eig$vectors %*% diag(exp(eig$values)) %*% t(eig$vectors)
}

#' An ellipsoid 6-vector from an unnormalized log-ellipsoid tensor.
#' 
#' The conversion is tuned so that the Frobenius inner product on matrices maps to the dot product on vectors.
#' @param logEll A real 3x3 matrix (symmetric).
#' @return A 6-dimensional real vector.
ellVectorFromLog <- function(logEll) {
  v1 <- sqrt(2) * logEll[1, 2]
  v2 <- sqrt(2) * logEll[1, 3]
  v3 <- sqrt(2) * logEll[2, 3]
  v4 <- logEll[1, 1]
  v5 <- logEll[2, 2]
  v6 <- logEll[3, 3]
  c(v1, v2, v3, v4, v5, v6)
}

#' An unnormalized log-ellipsoid tensor from an ellipsoid 6-vector.
#' 
#' The conversion is tuned so that the Frobenius inner product on matrices maps to the dot product on vectors.
#' @param vec A 6-dimensional real vector.
#' @return A real 3x3 matrix (symmetric).
ellLogFromVector <- function(vec) {
  l12 <- vec[[1]] / sqrt(2)
  l13 <- vec[[2]] / sqrt(2)
  l23 <- vec[[3]] / sqrt(2)
  l11 <- vec[[4]]
  l22 <- vec[[5]]
  l33 <- vec[[6]]
  matrix(c(l11, l12, l13, l12, l22, l23, l13, l23, l33), 3, 3)
}

#' An ellipsoid 5-vector from a normalized log-ellipsoid tensor.
#' 
#' The conversion is tuned so that the Frobenius inner product on matrices maps to the dot product on vectors.
#' @param logEll A real 3x3 matrix (symmetric, trace-zero).
#' @return A 5-dimensional real vector.
ellNormalizedVectorFromLog <- function(normLogEll) {
  v1 <- sqrt(2) * normLogEll[1, 2]
  v2 <- sqrt(2) * normLogEll[1, 3]
  v3 <- sqrt(2) * normLogEll[2, 3]
  v4 <- (normLogEll[2, 2] + normLogEll[1, 1]) * sqrt(1.5)
  v5 <- (normLogEll[2, 2] - normLogEll[1, 1]) / sqrt(2)
  c(v1, v2, v3, v4, v5)
}

#' An ellipsoid 5-vector from a normalized log-ellipsoid tensor.
#' 
#' The conversion is tuned so that the Frobenius inner product on matrices maps to the dot product on vectors.
#' @param logEll A real 3x3 matrix (symmetric, trace-zero).
#' @return A 5-dimensional real vector.
ellNormalizedVectorFromLogNew <- function(normLogEll) {
  v1 <- sqrt(2) * normLogEll[1, 2]
  v2 <- sqrt(2) * normLogEll[1, 3]
  v3 <- sqrt(2) * normLogEll[2, 3]
  v4 <- (normLogEll[1, 1] - normLogEll[2, 2]) / sqrt(2)
  v5 <- (normLogEll[1, 1] + normLogEll[2, 2] - 2 * normLogEll[3, 3]) / sqrt(6)
  c(v1, v2, v3, v4, v5)
}

#' A normalized log-ellipsoid tensor from an ellipsoid 5-vector.
#' 
#' The conversion is tuned so that the Frobenius inner product on matrices maps to the dot product on vectors.
#' @param vec A 5-dimensional real vector.
#' @return A real 3x3 matrix (symmetric, trace-zero).
ellLogFromNormalizedVectorNew <- function(vec) {
  l11 <- vec[[4]] / sqrt(2) + vec[[5]] / sqrt(6);
  l12 <- vec[[1]] / sqrt(2);
  l13 <- vec[[2]] / sqrt(2);
  l22 <- -vec[[4]] / sqrt(2) + vec[[5]] / sqrt(6);
  l23 <- vec[[3]] / sqrt(2);
  l33 <- -sqrt(2 / 3) * vec[[5]];
  matrix(c(l11, l12, l13, l12, l22, l23, l13, l23, l33), 3, 3)
}

#' A normalized log-ellipsoid tensor from an ellipsoid 5-vector.
#' 
#' The conversion is tuned so that the Frobenius inner product on matrices maps to the dot product on vectors.
#' @param vec A 5-dimensional real vector.
#' @return A real 3x3 matrix (symmetric, trace-zero).
ellLogFromNormalizedVector <- function(vec) {
  l11 <- vec[[4]] / sqrt(6) - vec[[5]] / sqrt(2)
  l22 <- vec[[4]] / sqrt(6) + vec[[5]] / sqrt(2)
  l33 <- -sqrt(2 / 3) * vec[[4]]
  l12 <- vec[[1]] / sqrt(2)
  l13 <- vec[[2]] / sqrt(2)
  l23 <- vec[[3]] / sqrt(2)
  matrix(c(l11, l12, l13, l12, l22, l23, l13, l23, l33), 3, 3)
}

#' An ellipsoid from an ellipsoid vector.
#' 
#' @param v An ellipsoid vector, either 5-dimensional (if normalized) or 6-dimensional (if not).
#' @return An ellipsoid.
ellEllipsoidFromVector <- function(v) {
  # Diagonalize the log ellipsoid tensor.
  if (length(v) == 5)
    logEll <- ellLogFromNormalizedVector(v)
  else
    logEll <- ellLogFromVector(v)
  eig <- eigen(logEll, symmetric=TRUE)
  # Everything else follows from that diagonalization.
  logA <- -0.5 * eig$values
  a <- exp(logA)
  rotation <- t(eig$vectors)
  if (det(rotation) < 0)
    rotation[3,] <- -rotation[3,]
  tensor <- eig$vectors %*% diag(exp(eig$values)) %*% t(eig$vectors)
  list(vector=v, a=a, logA=logA, rotation=rotation, tensor=tensor)
}

#' An ellipsoid from its orientation and the logarithms of its semi-axis lengths.
#' 
#' @param r A real 3x3 matrix (special orthogonal), with the semi-axis directions along its rows.
#' @param logA A real 3-dimensional vector. The logarithms of the semi-axis lengths, in order corresponding to the rows of r.
#' @return An ellipsoid.
ellEllipsoidFromRotationLogA <- function(r, logA, doNormalize=FALSE) {
  if (doNormalize)
    logA <- logA - sum(logA) / 3
  a <- exp(logA)
  tensor <- t(r) %*% diag(a^-2) %*% r
  logEll <- t(r) %*% diag(-2 * logA) %*% r
  if (doNormalize)
    v <- ellNormalizedVectorFromLog(logEll)
  else
    v <- ellVectorFromLog(logEll)
  list(vector=v, a=a, logA=logA, rotation=r, tensor=tensor)
}

ellTensorFromRotationA <- function(r, a) {
  t(r) %*% diag(a^-2) %*% r
}

ellRotationAFromTensor <- function(e) {
  eig <- eigen(e, symmetric=TRUE)
  a <- eig$values^(-1 / 2)
  r <- t(eig$vectors)
  if (det(r) < 0)
    r[3,] <- -r[3,]
  list(rotation=r, a=a)
}

ellEllipsoidFromTensor <- function(e, doNormalize=FALSE) {
  ra <- ellRotationAFromTensor(e)
  logA <- log(ra$a)
  if (doNormalize)
    logA <- logA - sum(logA) / 3
  logEll <- t(ra$rotation) %*% diag(-2 * logA) %*% ra$rotation
  if (doNormalize)
    v <- ellNormalizedVectorFromLog(logEll)
  else
    v <- ellVectorFromLog(logEll)
  list(vector=v, logA=logA, a=exp(logA), rotation=ra$rotation, tensor=e)
}

#' Ellipsoid orientation with axes in ascending or descending order.
#' 
#' @param rot A real 3x3 matrix (special orthogonal), with the semi-axis directions along its rows.
#' @param aOrLogA A real 3-dimensional vector. The semi-axis lengths or their logarithms, in order corresponding to the rows of rot.
#' @param descending A Boolean. If TRUE, then sort axes in descending order. If FALSE, then sort in ascending order.
#' @return A real 3x3 matrix (special orthogonal). This matrix equals the input matrix, with its rows reordered and possibly one row negated to maintain determinant 1.
ellAscendingRotation <- function(rot, aOrLogA, descending=FALSE) {
  ord <- order(aOrLogA, decreasing=descending)
  first <- rot[ord[[1]],]
  second <- rot[ord[[2]],]
  third <- rot[ord[[3]],]
  if (dot(cross(first, second), third) > 0)
    rbind(first, second, third)
  else
    rbind(first, second, -third)
}

ellDescendingRotation <- function(rot, aOrLogA) {
  ellAscendingRotation(rot, aOrLogA, descending=TRUE)
}

#' Ellipsoid orientation with axes ordered short, then long, then intermediate.
#' 
#' This function is useful for comparing ellipsoid orientations to foliation-lineation orientations.
#' @param rot A real 3x3 matrix (special orthogonal), with the semi-axis directions along its rows.
#' @param aOrLogA A real 3-dimensional vector. The semi-axis lengths or their logarithms, in order corresponding to the rows of rot.
#' @return A real 3x3 matrix (special orthogonal). This matrix equals the input matrix, with its rows reordered and possibly one row negated to maintain determinant 1.
ellPoleDirectionRotation <- function(rot, aOrLogA) {
  ord <- order(aOrLogA)
  first <- rot[ord[[1]],]
  second <- rot[ord[[3]],]
  third <- rot[ord[[2]],]
  if (dot(cross(first, second), third) > 0)
    rbind(first, second, third)
  else
    rbind(first, second, -third)
}



### SIZE AND SHAPE ###

#' Size-related tensor invariant of ellipsoids. Tantamount to volume.
#' 
#' The first invariant (trace) of log E^(-1 / 2), where E is the ellipsoid tensor. The volume of the ellipsoid is (4 pi / 3) * exp(size). The size is positive for large ellipsoids, zero for normalized ellipsoids, and negative for small ellipsoids.
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number (can achieve any real value).
ellSizeInvariant <- function(logs) {
  sum(logs)
}

#' Strain-related tensor invariant of ellipsoids. Tantamount to octahedral shear strain.
#'
#' The second invariant of log E^(-1 / 2), where E is the ellipsoid tensor. Equals zero for spheres. For normalized ellipsoids, this strain == -Es^2 / 2, where Es is the octahedral shear strain. 
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number <= 0.
ellStrainInvariant <- function(logs) {
  logs[[1]] * logs[[2]] + logs[[2]] * logs[[3]] + logs[[3]] * logs[[1]]
}

#' Shape-related tensor invariant of ellipsoids. An analogue of Lode's parameter.
#'
#' The third invariant (determinant) of log E^(-1 / 2), where E is the ellipsoid tensor. For normalized ellipsoids, this shape is positive for oblate ellipsoids, zero for 'plane strain' ellipsoids, and negative for prolate ellipsoids. In this sense it is analogous to (but not equal to or even tantamount to) Lode's parameter.
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number (can achieve any real value).
ellShapeInvariant <- function(logs) {
  logs[[1]] * logs[[2]] * logs[[3]]
}

#' Volume of an ellipsoid.
#' 
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number > 0.
ellVolume <- function(logs) {
  exp(sum(logs)) * 4 * pi / 3
}

#' Octahedral shear strain e_s.
#' 
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number >= 0.
ellOctahedralShearStrain <- function(logs) {
  mostOfIt <- (logs[[1]] - logs[[2]])^2 + (logs[[2]] - logs[[3]])^2 + (logs[[3]] - logs[[1]])^2
  sqrt(mostOfIt / 3)
}

# Lode's parameter nu.
#' 
#' nu is undefined for spheres, but we arbitrarily declare nu = 0 for them. Otherwise -1 <= nu <= 1. nu = -1 for prolate spheroids and nu = 1 for oblate spheroids.
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number (in the interval [-1, 1]), unless it fails.
ellLodeNu <- function(logs) {
  # Sort the logs so that l1 >= l2 >= l3.
  l1 <- max(logs)
  l3 <- min(logs)
  l2 <- sum(logs) - l1 - l3
  if (l1 == l3)
    0
  else
    (2 * l2 - l1 - l3) / (l1 - l3)
}

#' The statistic P_j of Jelinek (1981).
#' 
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number >= 1.
ellJelinekP <- function(logs) {
  # We used to have an '8' in the square root, where Jelinek had a '2', because we thought that he was working with eta_i = -2 l_i. But on 2016/03/10 Mike Jackson at the IRM told me that, according to all recent authors, AMS ellipsoids are 'magnitude ellipsoids'(e.g., Hrouda, 1982), whose semi-axis lengths are the principal susceptibilities, which are the eigenvalues of the susceptibility tensor. So it seems that eta_i = l_i now. And saying so reproduces the IRM's computed P_j.
  v <- logs - sum(logs) / 3
  exp(sqrt(2 * dot(v, v)))
}

#' Flinn's K measure of ellipsoid shape.
#' 
#' Fails in the case of a prolate spheroid or sphere. Zero in the case of an oblate spheroid (that is not a sphere).
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number >= 0, unless it fails.
ellFlinnK <- function(logs) {
  # Sort the logs so that l1 >= l2 >= l3.
  l1 <- max(logs)
  l3 <- min(logs)
  l2 <- sum(logs) - l1 - l3
  a1 <- exp(l1)
  a2 <- exp(l2)
  a3 <- exp(l3)
  (a1 / a2 - 1) / (a2 / a3 - 1)
}

#' The logs of an ellipsoid's semi-axis lengths, from three other measures of shape.
#' 
#' @param vEsNu A real 3D vector consisting of volume, octahedral shear strain, and Lode's parameter nu.
#' @return A real 3D vector, consisting of the logs of the ellipsoid's semi-axis lengths.
ellLogAFromVEsNu <- function(vEsNu) {
  # Invert the volume-logs relationship.
  sumOfLogs <- log(vEsNu[[1]] * 3 / (4 * pi))
  # logs1 = alpha + beta logs3.
  alpha <- sumOfLogs * 2 / (vEsNu[[3]] + 3)
  beta <- (vEsNu[[3]] - 3) / (vEsNu[[3]] + 3)
  # logs2 = gamma + delta logs3.
  gamma <- (1 - 2 / (vEsNu[[3]] + 3)) * sumOfLogs
  delta <- -1 - beta
  # Compute the coefficients of the quadratic.
  aa <- 2 - 2 * beta + 2 * beta^2 - 2 * delta - 2 * beta * delta + 2 * delta^2
  bb <- -2 * alpha + 4 * alpha * beta - 2 * gamma - 2 * beta * gamma - 2 * alpha * delta + 4 * gamma * delta
  cc <- 2 * alpha^2 - 2 * alpha * gamma + 2 * gamma^2 - 3 * vEsNu[[2]]^2
  # Solve the quadratic aa logs3^2 + bb logs3 + cc == 0 and back out the other logs.
  logs3 <- realQuadraticSolutions(aa, bb, cc)
  logs1 <- sapply(logs3, function(l3) {alpha + beta * l3})
  logs2 <- sapply(logs3, function(l3) {gamma + delta * l3})
  sols <- lapply(1:length(logs3), function(i) c(logs1[[i]], logs2[[i]], logs3[[i]]))
  # Choose the solution such that logs1 >= logs2 >= logs3, as required by nu.
  sols <- Filter(function(sol) {sol[[1]] >= sol[[2]] && sol[[2]] >= sol[[3]]}, sols)
  if (length(sols) != 1) {
    print("warning: ellLogsFromVEsNu: did not find one and only one solution as expected")
    print(sols)
  }
  sols[[1]]
}



### DESCRIPTIVE STATISTICS ###

#' Geometric mean of ellipsoids.
#' 
#' @param vectors A list of 5- or 6-dimensional real vectors.
#' @return An ellipsoid.
ellMean <- function(vectors) {
  ellEllipsoidFromVector(arithmeticMean(vectors))
}

#' Covariance matrix.
#' 
#' The vectors are automatically centered about their mean. The denominator is n - 1, not n.
#' @param vectors A list of 5- or 6-dimensional real vectors.
#' @return A 5x5 or 6x6 real matrix.
ellCovariance <- function(vectors) {
  var(t(simplify2array(vectors)))
}

#' Convenience shortcut to the square roots of the eigenvalues of the covariance matrix.
#' 
#' These 5 or 6 numbers quantify the dispersion of the given ellipsoids. Together they are the multivariate analogue of the sample standard deviation.
#' @param vectors A list of 5- or 6-dimensional real vectors.
#' @return A 5- or 6-dimensional vector, containing the eigenvalues of the variance.
ellCovarianceScalars <- function(vectors) {
  sqrt(eigen(ellCovariance(vectors), only.values=TRUE)$values)
}

#' Principal component analysis.
#' 
#' @param vectors A list of 5- or 6-dimensional real vectors.
#' @return See prcomp. For quick reference, as.numeric(t(pca$rotation) %*% (v - pca$center)) is the vector v translated and rotated to its corresponding vector of weights. That's probably what you want to plot.
ellPrincipalComponentAnalysis <- function(vectors) {
  prcomp(t(simplify2array(vectors)))
}



### INFERENCE ###

#' One-sample inference based on Hotelling's T2 test.
#' 
#' @param vectors A list of 5- or 6-dimensional real vectors.
#' @param hypoth A 5- or 6-dimensional real vector, respectively.
#' @param fOrChi Character. Should be 'f' or 'chi'.
#' @return See HotellingsT2.
ellHotellingT2Inference <- function(vectors, hypoth, fOrChi="f") {
  HotellingsT2(X=t(simplify2array(vectors)), mu=hypoth, test=fOrChi)
}

# 32 or 64 corners of a crude confidence region.
ellCICombinationVectors <- function(ci) {
  if (!class(ci) == "matrix")
    list(ci[[1]], ci[[2]])
  else {
    recursive <- ellCICombinationVectors(ci[,(2:ncol(ci))])
    c(lapply(recursive, function(rec) c(ci[[1, 1]], rec)),
      lapply(recursive, function(rec) c(ci[[2, 1]], rec)))
  }
}

# Generates approximately 7^(dim - 1) points on the (dim - 1)-dimensional unit sphere in dim-dimensional Euclidean space.
ellHighSphereVectors <- function(ambientDimension, numSamples=7) {
  if (ambientDimension == 0)
    list()
  else if (ambientDimension == 1)
    list(1, -1)
  else if (ambientDimension == 2)
    lapply(0:(numSamples - 1),
           function(i) {a <- (i * 2 * pi + 1) / numSamples; c(sin(a), cos(a))})
  else {
    recursive <- ellHighSphereVectors(ambientDimension - 1, numSamples)
    unlist(lapply(
      recursive,
      function(rec) lapply(0:(numSamples - 1),
                           function(i) {a <- (i * 2 * pi + 1) / numSamples; c(sin(a) * rec, cos(a))})),
      recursive=FALSE, use.names=FALSE)
  }
}

#' A sampling of points on an ellipsoid in high-dimensional space.
#' 
#' d is the dimension of the ambient space. The arguments for this function typically come out of ellBootstrapInference, where d == 5 or d == 6.
#' @param covarInv A real dxd matrix (symmetric, positive-definite).
#' @param center A real dD vector.
#' @param level A real number. Typically q095^2 from ellBootstrapInference.
#' @param numSamples A real number (positive integer). Roughly, the number of samples per dimension.
#' @return A list of dD real vectors.
ellHighEllipsoidVectors <- function(covarInv, center, level, numSamples=7) {
  eig <- eigen(covarInv, symmetric=TRUE)
  q <- eig$vectors
  a <- sqrt(level) * eig$values^(-0.5)
  sphere <- ellHighSphereVectors(ncol(covarInv), numSamples)
  lapply(sphere, function(v) as.numeric((q %*% (a * v)) + center))
}

ellMahalanobisDistance <- function(ss, sBar) {
  vs <- lapply(ss, function(s) {s - sBar})
  covar <- arithmeticMean(lapply(vs, function(v) {outer(v, v)}))
  covarInv <- solve(covar)
  f <- function(s) {
    v <- s - sBar
    sqrt(v %*% covarInv %*% v)
  }
  f
}

# Confidence region based on percentiles of Mahalanobis distance.
ellMahalanobisPercentiles <- function(ss, sBar) {
  vs <- lapply(ss, function(s) {s - sBar})
  covar <- arithmeticMean(lapply(vs, function(v) {outer(v, v)}))
  covarInv <- solve(covar)
  norms <- sapply(vs, function(v) {sqrt(v %*% covarInv %*% v)})
  empiricalCDF <- ecdf(norms)
  # Build the p-value function.
  f <- function(s) {
    v <- s - sBar
    1 - empiricalCDF(sqrt(v %*% covarInv %*% v))
  }
  # Compute a few popular percentiles.
  qs <- quantile(norms, probs=c(0.00, 0.25, 0.50, 0.75, 0.95, 0.99, 1.00), names=FALSE)
  list(pvalue=f, center=sBar, covarInv=covarInv,
       q000=qs[[1]], q025=qs[[2]], q050=qs[[3]], q075=qs[[4]], q095=qs[[5]], q099=qs[[6]], q100=qs[[7]])
}

#' One-sample inference based on bootstrapping.
#' 
#' @param vectors A list of 5- or 6-dimensional real vectors.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @return A list with members $pvalue, $center, $covarInv, $bootstraps, $q000, $q025, $q050, $q075, $q095, $q099, $q100. $pvalue is a function with input a 5D or 6D vector and output a real number, the p-value for the null hypothesis that the mean is that vector. $center is the mean of $bootstraps, which are the bootstraps. $covarInv is their inverse covariance matrix at the mean. The $qxxx values are percentiles of Mahalanobis distance among the bootstraps.
ellBootstrapInference <- function(vectors, numBoots=1000) {
  boots <- replicate(numBoots, arithmeticMean(sample(vectors, length(vectors), replace=TRUE)), simplify=FALSE)
  bootMean <- arithmeticMean(boots)
  infer <- ellMahalanobisPercentiles(boots, bootMean)
  infer$bootstraps <- boots
  infer
}

#' Two-sample inference based on Hotelling's T2 test.
#' 
#' @param firsts A list of 5- or 6-dimensional real vectors.
#' @param seconds A list of 5- or 6-dimensional real vectors.
#' @param fOrChi Character. Should be 'f' or 'chi'.
#' @return See HotellingsT2.
ellTwoSampleHotellingT2Inference <- function(firsts, seconds, fOrChi="f") {
  HotellingsT2(X=t(simplify2array(firsts)), Y=t(simplify2array(seconds)), test=fOrChi)
}

#' Two-sample inference based on bootstrapping.
#' 
#' @param firsts A list of 5- or 6-dimensional real vectors.
#' @param seconds A list of 5- or 6-dimensional real vectors.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @return A list with members $pvalue, $center, $covarInv, $bootstraps, $q000, $q025, $q050, $q075, $q095, $q099, $q100. $pvalue is a function with input a 5D or 6D vector and output a real number, the p-value for the null hypothesis that second-mean - first-mean is that vector. $center is the mean of $bootstraps, which are the bootstraps. $covarInv is their inverse covariance matrix at the mean. The $qxxx values are percentiles of Mahalanobis distance among the bootstraps.
ellTwoSampleBootstrapInference <- function(firsts, seconds, numBoots=1000) {
  f <- function() {
    firstMean <- arithmeticMean(sample(firsts, length(firsts), replace=TRUE))
    secondMean <- arithmeticMean(sample(seconds, length(seconds), replace=TRUE))
    secondMean - firstMean
  }
  boots <- replicate(numBoots, f(), simplify=FALSE)
  bootMean <- arithmeticMean(boots)
  infer <- ellMahalanobisPercentiles(boots, bootMean)
  infer$bootstraps <- boots
  infer
}



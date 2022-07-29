


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### WELLNER (1979) TWO-SAMPLE TEST ###

#' Wellner's T-statistic (Wellner, 1979, Example 1b), which measures how different two sets of lines are.
#' 
#' @param xs A list of lines.
#' @param ys A list of lines.
#' @return A real number (nonnegative). Zero if the two data sets are identical.
lineWellner <- function(xs, ys) {
  # Prepare the number of xs, the number of ys, and the dimension of the hypersphere.
  m <- length(xs)
  n <- length(ys)
  p <- length(xs[[1]]) - 1
  # Compute m Tx and n Ty.
  f <- function(v) {outer(v, v)}
  mTx <- Reduce("+", lapply(xs, f))
  nTy <- Reduce("+", lapply(ys, f))
  # Compute Wellner's statistic.
  txMinusTy <- mTx / m - nTy / n
  tr(txMinusTy %*% txMinusTy) * m * n * (p + 1) * (p + 3) / (2 * (m + n))
}

# Also computes Wellner's T-statistic, with inputs presented in a manner suitable for fast permutation testing.
lineWellnerShortcut <- function(xys, choices, mTxPlusnTy) {
  # Prepare the number of xs, the number of ys, and the dimension of the hypersphere.
  m <- length(choices)
  n <- length(xys) - m
  p <- length(xys[[1]]) - 1
  # Compute m Tx and n Ty.
  mTx <- Reduce("+", lapply(choices, function(i) {outer(xys[[i]], xys[[i]])}))
  nTy <- mTxPlusnTy - mTx
  # Compute Wellner's statistic.
  txMinusTy <- mTx / m - nTy / n
  tr(txMinusTy %*% txMinusTy) * m * n * (p + 1) * (p + 3) / (2 * (m + n))
}

#' Two-sample test, based on permutations and Wellner's T-statistic (Wellner, 1979).
#' 
#' @param xs A list of lines.
#' @param ys A list of lines.
#' @param numPerms A real number (positive integer). The number of permutations, say 1,000 or 10,000.
#' @return A real number, between 0 and 1 inclusive. The fraction of tests in which T exceeds the original T for the data. You can interpret this as a p-value for the null hypothesis that the two populations are identical (not just that their means are identical). In other words, small values of p indicate that the distinction between the two populations is meaningful.
lineWellnerInference <- function(xs, ys, numPerms) {
  # Precompute.
  m <- length(xs)
  n <- length(ys)
  xys <- c(xs, ys)
  mTxPlusnTy <- Reduce("+", lapply(xys, function(v) {outer(v, v)}))
  # Compute Wellner's statistic for the actual data.
  t <- lineWellnerShortcut(xys, 1:m, mTxPlusnTy)
  # Compute Wellner's statistic for permuted data, numPerms times.
  ts <- replicate(numPerms, lineWellnerShortcut(xys, sample.int(m + n, m), mTxPlusnTy))
  # What proportion of permuted results are more extreme than the one observed?
  greaterThan <- Filter(function (u) {u > t}, ts)
  length(greaterThan) / numPerms
}

# Tests.
lineWellnerInferenceExperimentFuncAs <- function(eastValues, westValues, sameCenter, sameDirection) {
  function() {
    # Choose the centers. If sameDirection, then force sameCenter.
    east3 <- rayUniform()
    if (sameCenter || sameDirection)
      west3 <- east3
    else
      west3 <- rayUniform()
    # Choose the directions.
    east1 <- rayOrthogonalUniform(east3)
    if (sameDirection)
      west1 <- east1
    else
      west1 <- rayOrthogonalUniform(west3)
    eastRot <- cbind(east1, cross(east3, east1), east3)
    westRot <- cbind(west1, cross(west3, west1), west3)
    eastA <- eastRot %*% diag(eastValues) %*% t(eastRot)
    westA <- westRot %*% diag(westValues) %*% t(westRot)
    list(eastA, westA)
  }
}
lineWellnerInferenceExperiment <- function(eastValues, westValues, sameCenter, sameDirection, eastN, westN, numPerms, numTrials) {
  f <- function(i) {
    print(i / numTrials)
    # Choose the centers. If sameDirection, then force sameCenter.
    east3 <- rayUniform()
    if (sameCenter || sameDirection)
      west3 <- east3
    else
      west3 <- rayUniform()
    # Choose the directions.
    east1 <- rayOrthogonalUniform(east3)
    if (sameDirection)
      west1 <- east1
    else
      west1 <- rayOrthogonalUniform(west3)
    eastRot <- cbind(east1, cross(east3, east1), east3)
    westRot <- cbind(west1, cross(west3, west1), west3)
    eastA <- eastRot %*% diag(eastValues) %*% t(eastRot)
    westA <- westRot %*% diag(westValues) %*% t(westRot)
    easts <- lineBingham(eastA, eastN)
    wests <- lineBingham(westA, westN)
    lineWellnerInference(easts, wests, numPerms=numPerms)
  }
  ps <- sapply(1:numTrials, f)
  pHat <- sum(ps > 0.05) / numTrials
  c(pHat, 2 * standardErrorProportion(numTrials, pHat))
}
#lineWellnerInferenceExperiment(
#  eastValues=c(1, 0, -1), westValues=c(1, 0, -1), sameCenter=TRUE, 
#  sameDirection=TRUE, eastN=10, westN=10, numPerms=1000, numTrials=2000)
# 0.94550000 0.01015182


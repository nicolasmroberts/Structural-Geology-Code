


# Copyright 2018 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### WELLNER (1979) TWO-SAMPLE TEST ###

#' Wellner's Rayleigh-style T-statistic (Wellner, 1979, Example 1a), which measures how different two sets of rays are.
#' 
#' @param xs A list of rays.
#' @param ys A list of rays.
#' @return A real number (nonnegative). Zero if the two data sets are identical.
rayWellner <- function(xs, ys) {
  # Prepare the number of xs, the number of ys, and the dimension of the hypersphere.
  m <- length(xs)
  n <- length(ys)
  p <- length(xs[[1]]) - 1
  # Compute Wellner's statistic.
  diff <- arithmeticMean(xs) - arithmeticMean(ys)
  dot(diff, diff) * (p + 1) * m * n / (m + n)
}

# Also computes Wellner's T-statistic, with inputs presented in a manner suitable for fast permutation testing.
rayWellnerShortcut <- function(xys, choices, rxPlusRy) {
  # Prepare the number of xs and ys and the dimension of the hypersphere.
  m <- length(choices)
  n <- length(xys) - m
  p <- length(xys[[1]]) - 1
  # Compute Rx and Ry.
  rx <- Reduce("+", lapply(choices, function(i) xys[[i]]))
  ry <- rxPlusRy - rx
  # Compute Wellner's statistic.
  diff <- rx / m - ry / n
  dot(diff, diff) * (p + 1) * m * n / (m + n)
}

#' Two-sample test, based on permutations and Wellner's Rayleigh-style T-statistic (Wellner, 1979, Example 1a).
#' 
#' Assumes large sample sizes, specifically that choose(nx + ny, nx) >> numPerms. For small sample sizes, see rayWellnerExactInference.
#' @param xs A list of rays.
#' @param ys A list of rays.
#' @param numPerms A real number (positive integer). The number of permutations, say 1,000 or 10,000.
#' @return A real number, between 0 and 1 inclusive. The fraction of tests in which T exceeds the original T for the data. You can interpret this as a p-value for the null hypothesis that the two populations are identical (not just that their means are identical). In other words, small values of p indicate that the distinction between the two populations is meaningful.
rayWellnerInference <- function(xs, ys, numPerms) {
  # Precompute.
  m <- length(xs)
  n <- length(ys)
  xys <- c(xs, ys)
  rxPlusRy <- Reduce("+", xys)
  # Compute Wellner's statistic for the actual data.
  t <- rayWellnerShortcut(xys, 1:m, rxPlusRy)
  # Compute Wellner's statistic for permuted data, numPerms times.
  ts <- replicate(numPerms, rayWellnerShortcut(xys, sample.int(m + n, m), rxPlusRy))
  # What proportion of permuted results are more extreme than the one observed?
  greaterThan <- Filter(function (u) {u > t}, ts)
  length(greaterThan) / numPerms
}

#' Two-sample test, based on permutations and Wellner's Rayleigh-style T-statistic (Wellner, 1979, Example 1a).
#' 
#' Assumes small sample sizes. Deterministically generates all choose(nx + ny, nx) reassignments of the data to the two groups. For large sample sizes, see rayWellnerInference.
#' @param xs A list of rays.
#' @param ys A list of rays.
#' @return A real number, between 0 and 1 inclusive. The fraction of tests in which T exceeds the original T for the data. You can interpret this as a p-value for the null hypothesis that the two populations are identical (not just that their means are identical). In other words, small values of p indicate that the distinction between the two populations is meaningful.
rayWellnerExactInference <- function(xs, ys) {
  # Precompute.
  m <- length(xs)
  n <- length(ys)
  xys <- c(xs, ys)
  rxPlusRy <- Reduce("+", xys)
  # Compute Wellner's statistic for the actual data.
  t <- rayWellnerShortcut(xys, 1:m, rxPlusRy)
  # Compute Wellner's statistic for all possible data permutations.
  rows <- combs(1:(m + n), m)
  ts <- sapply(
    1:nrow(rows), 
    function(i) rayWellnerShortcut(xys, rows[i,], rxPlusRy))
  # What proportion of permuted results are more extreme than the one observed?
  greaterThan <- Filter(function (u) {u > t}, ts)
  length(greaterThan) / choose(m + n, m)
}

# Experiments on coverage rates with equal means.
rayWellnerInferenceExperiment <- function(n1, kappa1, n2, kappa2, numPerms, numTrials) {
  if (numPerms < choose(n1 + n2, n1)) {
    f <- function(i) {
      print(i / numTrials)
      mu <- rayUniform()
      data1 <- rayFisher(mu, kappa1, n1)
      data2 <- rayFisher(mu, kappa2, n2)
      rayWellnerInference(data1, data2, numPerms)
    }
  } else {
    f <- function(i) {
      print(i / numTrials)
      mu <- rayUniform()
      data1 <- rayFisher(mu, kappa1, n1)
      data2 <- rayFisher(mu, kappa2, n2)
      rayWellnerExactInference(data1, data2)
    }
  }
  ps <- sapply(1:numTrials, f)
  pHat <- sum(ps > 0.05) / numTrials
  c(pHat, 2 * standardErrorProportion(numTrials, pHat))
}

# Experiments with equal means, kappa1 = kappa2 = 1, varying sample size. The 
# test should accept (fail to reject) 95% of the time, and indeed it does.
#rayWellnerInferenceExperiment(3, 1, 3, 1, 1000, 2000)#0.90500000 0.01311297
#rayWellnerInferenceExperiment(3, 1, 10, 1, 1000, 2000)#0.94500000 0.01019559
#rayWellnerInferenceExperiment(3, 1, 30, 1, 1000, 2000)#0.953500000 0.009416767
#rayWellnerInferenceExperiment(3, 1, 100, 1, 1000, 2000)#0.93750000 0.01082532
#rayWellnerInferenceExperiment(10, 1, 10, 1, 1000, 2000)#0.957500000 0.009021502
#rayWellnerInferenceExperiment(10, 1, 30, 1, 1000, 2000)#0.952500000 0.009512492
#rayWellnerInferenceExperiment(10, 1, 100, 1, 1000, 2000)#0.955000000 0.009270922
#rayWellnerInferenceExperiment(30, 1, 30, 1, 1000, 2000)#0.94450000 0.01023912
#rayWellnerInferenceExperiment(30, 1, 100, 1, 1000, 2000)#0.94450000 0.01023912
#rayWellnerInferenceExperiment(100, 1, 100, 1, 1000, 2000)#0.955000000 0.009270922

# Experiments with equal means, kappa1 = kappa2 = 10, varying sample size.
#rayWellnerInferenceExperiment(3, 10, 3, 10, 1000, 2000)#0.89350000 0.01379549
#rayWellnerInferenceExperiment(3, 10, 10, 10, 1000, 2000)#0.93200000 0.01125842
#rayWellnerInferenceExperiment(3, 10, 30, 10, 1000, 2000)#0.95050000 0.00970049
#rayWellnerInferenceExperiment(3, 10, 100, 10, 1000, 2000)#0.94400000 0.01028241
#rayWellnerInferenceExperiment(10, 10, 10, 10, 1000, 2000)#0.94800000 0.00992935
#rayWellnerInferenceExperiment(10, 10, 30, 10, 1000, 2000)#0.94400000 0.01028241
#rayWellnerInferenceExperiment(10, 10, 100, 10, 1000, 2000)#0.953000000 0.009464777
#rayWellnerInferenceExperiment(30, 10, 30, 10, 1000, 2000)#0.949500000 0.009792829
#rayWellnerInferenceExperiment(30, 10, 100, 10, 1000, 2000)#0.93450000 0.01106433
#rayWellnerInferenceExperiment(100, 10, 100, 10, 1000, 2000)#0.94450000 0.01023912

# Experiments with equal means, kappa1 = kappa2 = 100, varying sample size.
#rayWellnerInferenceExperiment(3, 100, 3, 100, 1000, 2000)#0.9380000 0.0107848
#rayWellnerInferenceExperiment(3, 100, 10, 100, 1000, 2000)#0.953000000 0.009464777
#rayWellnerInferenceExperiment(3, 100, 30, 100, 1000, 2000)#0.94450000 0.01023912
#rayWellnerInferenceExperiment(3, 100, 100, 100, 1000, 2000)#0.955000000 0.009270922
#rayWellnerInferenceExperiment(10, 100, 10, 100, 1000, 2000)#0.954000000 0.009368458
#rayWellnerInferenceExperiment(10, 100, 30, 100, 1000, 2000)#0.959500000 0.008815866
#rayWellnerInferenceExperiment(10, 100, 100, 100, 1000, 2000)#0.955500000 0.009221686
#rayWellnerInferenceExperiment(30, 100, 30, 100, 1000, 2000)#0.951500000 0.009607055
#rayWellnerInferenceExperiment(30, 100, 100, 100, 1000, 2000)#0.94800000 0.00992935
#rayWellnerInferenceExperiment(100, 100, 100, 100, 1000, 2000)#0.94300000 0.01036832

# Experiments with equal means and sample sizes but differing kappas. We'd like 
# the test to reject, but that's not guaranteed (even asymptotically) by 
# Wellner's consistency result. Nevertheless, it does pretty well, except for 
# small sample sizes and high concentrations. Then the test is a bit too 
# conservative (fails to reject too often).
#rayWellnerInferenceExperiment(3, 1, 3, 10, 1000, 2000)#0.71150000 0.02026168
#rayWellnerInferenceExperiment(3, 1, 3, 100, 1000, 2000)#0.52750000 0.02232683
#rayWellnerInferenceExperiment(3, 10, 3, 100, 1000, 2000)#0.85400000 0.01579139
#rayWellnerInferenceExperiment(10, 1, 10, 10, 1000, 2000)#0.2540000 0.0194671
#rayWellnerInferenceExperiment(10, 1, 10, 100, 1000, 2000)#0.028500000 0.007441472
#rayWellnerInferenceExperiment(10, 10, 10, 100, 1000, 2000)#0.91300000 0.01260405
#rayWellnerInferenceExperiment(30, 1, 30, 10, 1000, 2000)#0 0
#rayWellnerInferenceExperiment(30, 1, 30, 100, 1000, 2000)#0 0
#rayWellnerInferenceExperiment(30, 10, 30, 100, 1000, 2000)#0.83450000 0.01661985
#rayWellnerInferenceExperiment(100, 1, 100, 10, 1000, 2000)#0 0
#rayWellnerInferenceExperiment(100, 1, 100, 100, 1000, 2000)#0 0
#rayWellnerInferenceExperiment(100, 10, 100, 100, 1000, 2000)#0.005500000 0.003307491

# Experiments on coverage rates with unequal means. Wellner's consistency 
# result says that the test should reject the null hypothesis on these almost 
# always (as the sample size goes to infinity).
rayWellnerInferenceUnequalExperiment <- function(n1, kappa1, n2, kappa2, numPerms, numTrials) {
  if (numPerms < choose(n1 + n2, n1)) {
    f <- function(i) {
      print(i / numTrials)
      data1 <- rayFisher(rayUniform(), kappa1, n1)
      data2 <- rayFisher(rayUniform(), kappa2, n2)
      rayWellnerInference(data1, data2, numPerms)
    }
  } else {
    f <- function(i) {
      print(i / numTrials)
      data1 <- rayFisher(rayUniform(), kappa1, n1)
      data2 <- rayFisher(rayUniform(), kappa2, n2)
      rayWellnerExactInference(data1, data2)
    }
  }
  ps <- sapply(1:numTrials, f)
  pHat <- sum(ps > 0.05) / numTrials
  c(pHat, 2 * standardErrorProportion(numTrials, pHat))
}

# Experiments with unequal means, kappa1 = kappa2 = 1, varying sample size. The 
# test should reject most of the time.
#rayWellnerInferenceUnequalExperiment(3, 1, 3, 1, 1000, 2000)#0.84650000 0.01612065
#rayWellnerInferenceUnequalExperiment(3, 1, 10, 1, 1000, 2000)#0.83600000 0.01655923
#rayWellnerInferenceUnequalExperiment(3, 1, 30, 1, 1000, 2000)#0.78250000 0.01844959
#rayWellnerInferenceUnequalExperiment(3, 1, 100, 1, 1000, 2000)#0.78300000 0.01843426
#rayWellnerInferenceUnequalExperiment(10, 1, 10, 1, 1000, 2000)#0.71500000 0.02018787
#rayWellnerInferenceUnequalExperiment(10, 1, 30, 1, 1000, 2000)#0.58800000 0.02201163
#rayWellnerInferenceUnequalExperiment(10, 1, 100, 1, 1000, 2000)#0.50900000 0.02235706
#rayWellnerInferenceUnequalExperiment(30, 1, 30, 1, 1000, 2000)#0.36350000 0.02151129
#rayWellnerInferenceUnequalExperiment(30, 1, 100, 1, 1000, 2000)#0.2350000 0.0189618
#rayWellnerInferenceUnequalExperiment(100, 1, 100, 1, 1000, 2000)#0.10700000 0.01382396

# Experiments with unequal means, kappa1 = kappa2 = 10, varying sample size.
#rayWellnerInferenceUnequalExperiment(3, 10, 3, 10, 1000, 2000)#0.11650000 0.01434767
#rayWellnerInferenceUnequalExperiment(3, 10, 10, 10, 1000, 2000)#0.07700000 0.01192233
#rayWellnerInferenceUnequalExperiment(3, 10, 30, 10, 1000, 2000)#0.06400000 0.01094568
#rayWellnerInferenceUnequalExperiment(3, 10, 100, 10, 1000, 2000)#0.05200000 0.00992935
#rayWellnerInferenceUnequalExperiment(10, 10, 10, 10, 1000, 2000)#0.039000000 0.008657829
#rayWellnerInferenceUnequalExperiment(10, 10, 30, 10, 1000, 2000)#0.02000000 0.00626099
#rayWellnerInferenceUnequalExperiment(10, 10, 100, 10, 1000, 2000)#0.019000000 0.006105571
#rayWellnerInferenceUnequalExperiment(30, 10, 30, 10, 1000, 2000)#0.011500000 0.004768176
#rayWellnerInferenceUnequalExperiment(30, 10, 100, 10, 1000, 2000)#0.007500000 0.003858432
#rayWellnerInferenceUnequalExperiment(100, 10, 100, 10, 1000, 2000)#0.002000000 0.001997999

# Experiments with unequal means, kappa1 = kappa2 = 100, varying sample size.
#rayWellnerInferenceUnequalExperiment(3, 100, 3, 100, 1000, 2000)#0.011500000 0.004768176
#rayWellnerInferenceUnequalExperiment(3, 100, 10, 100, 1000, 2000)#0.010000000 0.004449719
#rayWellnerInferenceUnequalExperiment(3, 100, 30, 100, 1000, 2000)#0.004000000 0.002822765
#rayWellnerInferenceUnequalExperiment(3, 100, 100, 100, 1000, 2000)#0.004000000 0.002822765
#rayWellnerInferenceUnequalExperiment(10, 100, 10, 100, 1000, 2000)#0.003500000 0.002641117
#rayWellnerInferenceUnequalExperiment(10, 100, 30, 100, 1000, 2000)#0.002500000 0.002233271
#rayWellnerInferenceUnequalExperiment(10, 100, 100, 100, 1000, 2000)#0.002500000 0.002233271
#rayWellnerInferenceUnequalExperiment(30, 100, 30, 100, 1000, 2000)#0.002000000 0.001997999
#rayWellnerInferenceUnequalExperiment(30, 100, 100, 100, 1000, 2000)#0.001500000 0.001730751
#rayWellnerInferenceUnequalExperiment(100, 100, 100, 100, 1000, 2000)#0.001000000 0.001413506

# Experiments with unequal means and unequal kappas.
#rayWellnerInferenceUnequalExperiment(3, 1, 3, 10, 1000, 2000)#0.4380000 0.0221881
#rayWellnerInferenceUnequalExperiment(3, 1, 3, 100, 1000, 2000)#0.25600000 0.01951738
#rayWellnerInferenceUnequalExperiment(3, 10, 3, 100, 1000, 2000)#0.05450000 0.01015182
#rayWellnerInferenceUnequalExperiment(10, 1, 10, 10, 1000, 2000)#0.032500000 0.007930164
#rayWellnerInferenceUnequalExperiment(10, 1, 10, 100, 1000, 2000)#0.002500000 0.002233271
#rayWellnerInferenceUnequalExperiment(10, 10, 10, 100, 1000, 2000)#0.021500000 0.006486563
#rayWellnerInferenceUnequalExperiment(30, 1, 30, 10, 1000, 2000)#0 0
#rayWellnerInferenceUnequalExperiment(30, 1, 30, 100, 1000, 2000)#0 0
#rayWellnerInferenceUnequalExperiment(30, 10, 30, 100, 1000, 2000)#0.003000000 0.002445813
#rayWellnerInferenceUnequalExperiment(100, 1, 100, 10, 1000, 2000)#0 0
#rayWellnerInferenceUnequalExperiment(100, 1, 100, 100, 1000, 2000)#0 0
#rayWellnerInferenceUnequalExperiment(100, 10, 100, 100, 1000, 2000)#0 0






# Copyright 2018 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



### WATSON (1983) TWO-SAMPLE TEST ###

#' The V_w statistic of Watson (1983, Eq. 4.4).
#' 
#' See Tauxe (2010, Section 11.3.4 and Appendix C.2.1). This statistic measures how different the sample means of two data sets are. (It is always non-negative. It is zero when the two sample means are equal.) It exists to facilitate the two-sample test rayFisherVwInference.
#' @param r1s A list of rays.
#' @param r2s A list of rays.
#' @return A real number. The statistic V_w.
rayFisherVw <- function(r1s, r2s) {
  # Compute some statistics for the first data set.
  bigX1Bar <- Reduce("+", r1s)
  bigR1 <- sqrt(dot(bigX1Bar, bigX1Bar))
  n1 <- length(r1s)
  k1 <- (n1 - 1) / (n1 - bigR1)
  # Compute the same statistics for the second data set.
  bigX2Bar <- Reduce("+", r2s)
  bigR2 <- sqrt(dot(bigX2Bar, bigX2Bar))
  n2 <- length(r2s)
  k2 <- (n2 - 1) / (n2 - bigR2)
  # Combine.
  bigXHat <- k1 * bigX1Bar + k2 * bigX2Bar
  bigRw <- sqrt(dot(bigXHat, bigXHat))
  bigSw <- k1 * bigR1 + k2 * bigR2
  bigVw <- 2 * (bigSw - bigRw)
  bigVw
}

#' Two-sample test based on V_w-statistic of Watson (1983).
#' 
#' The idea is to bootstrap rayFisherVw (parametrically, from the Fisher distribution) and compare the original V_w to the percentiles of the bootstrapped V_w.
#' @param r1s A list of rays.
#' @param r2s A list of rays.
#' @param numBoots A real number (positive integer). 
#' @return A p-value for the null hypothesis that the two populations have the same mean.
rayFisherVwInference <- function(r1s, r2s, numBoots) {
  # Compute some statistics for the first data set.
  bigX1Bar <- Reduce("+", r1s)
  bigR1 <- sqrt(dot(bigX1Bar, bigX1Bar))
  n1 <- length(r1s)
  k1 <- (n1 - 1) / (n1 - bigR1)
  # Compute the same statistics for the second data set.
  bigX2Bar <- Reduce("+", r2s)
  bigR2 <- sqrt(dot(bigX2Bar, bigX2Bar))
  n2 <- length(r2s)
  k2 <- (n2 - 1) / (n2 - bigR2)
  # Compute the sample V_w and the bootstrapped V_w.
  vw <- rayFisherVw(r1s, r2s)
  mu <- rayUniform()
  vws <- replicate(
    numBoots, 
    rayFisherVw(rayFisher(mu, k1, n1), rayFisher(mu, k2, n2)))
  empirical <- ecdf(vws)
  1 - empirical(vw)
}

# Experiments on coverage rates using equal means.
rayFisherVwInferenceExperiment <- function(n1, kappa1, n2, kappa2, numBoots, numTrials) {
  f <- function(i) {
    if (i %% 250 == 0)
      print(date())
    mu <- rayUniform()
    data1 <- rayFisher(mu, kappa1, n1)
    data2 <- rayFisher(mu, kappa2, n2)
    rayFisherVwInference(data1, data2, numBoots)
  }
  ps <- sapply(1:numTrials, f)
  pHat <- sum(ps > 0.05) / numTrials
  c(pHat, 2 * standardErrorProportion(numTrials, pHat))
}

# Experiments with equal means, kappa1 = kappa2 = 1, varying sample size. The 
# test should accept (fail to reject) 95% of the time.
#rayFisherVwInferenceExperiment(3, 1, 3, 1, 1000, 2000)#0.956000000 0.009172132
#rayFisherVwInferenceExperiment(3, 1, 10, 1, 1000, 2000)#0.92350000 0.01188678
#rayFisherVwInferenceExperiment(3, 1, 30, 1, 1000, 1000)#0.90400000 0.01863159
#rayFisherVwInferenceExperiment(3, 1, 100, 1, 1000, 1000)#0.91100000 0.01800878
#rayFisherVwInferenceExperiment(10, 1, 10, 1, 1000, 1000)#0.90000000 0.01897367
#rayFisherVwInferenceExperiment(10, 1, 30, 1, 1000, 1000)#0.91000000 0.01809972
#rayFisherVwInferenceExperiment(10, 1, 100, 1, 1000, 1000)#0.90600000 0.01845687
#rayFisherVwInferenceExperiment(30, 1, 30, 1, 1000, 1000)#0.91200000 0.01791714
#rayFisherVwInferenceExperiment(30, 1, 100, 1, 1000, 1000)#0.91400000 0.01773178
#rayFisherVwInferenceExperiment(100, 1, 100, 1, 1000, 1000)#0.91600000 0.01754355

# Experiments with equal means, kappa1 = kappa2 = 10, varying sample size.
#rayFisherVwInferenceExperiment(3, 10, 3, 10, 1000, 1000)#0.95900000 0.01254097
#rayFisherVwInferenceExperiment(3, 10, 10, 10, 1000, 1000)#0.93500000 0.01559166
#rayFisherVwInferenceExperiment(3, 10, 30, 10, 1000, 1000)#0.93000000 0.01613691
#rayFisherVwInferenceExperiment(3, 10, 100, 10, 1000, 1000)#0.92900000 0.01624303
#rayFisherVwInferenceExperiment(10, 10, 10, 10, 1000, 1000)#0.95000000 0.01378405
#rayFisherVwInferenceExperiment(10, 10, 30, 10, 1000, 1000)#0.93600000 0.01547953
#rayFisherVwInferenceExperiment(10, 10, 100, 10, 1000, 1000)#0.94600000 0.01429461
#rayFisherVwInferenceExperiment(30, 10, 30, 10, 1000, 1000)#0.94800000 0.01404222
#rayFisherVwInferenceExperiment(30, 10, 100, 10, 1000, 1000)#0.95300000 0.01338522
#rayFisherVwInferenceExperiment(100, 10, 100, 10, 1000, 1000)#0.94300000 0.01466301

# Experiments with equal means, kappa1 = kappa2 = 100, varying sample size.
#rayFisherVwInferenceExperiment(3, 100, 3, 100, 1000, 1000)#0.96000000 0.01239355
#rayFisherVwInferenceExperiment(3, 100, 10, 100, 1000, 1000)#0.92700000 0.01645248
#rayFisherVwInferenceExperiment(3, 100, 30, 100, 1000, 1000)#0.93100000 0.01602985
#rayFisherVwInferenceExperiment(3, 100, 100, 100, 1000, 1000)#0.95100000 0.01365269
#rayFisherVwInferenceExperiment(10, 100, 10, 100, 1000, 1000)#0.95900000 0.01254097
#rayFisherVwInferenceExperiment(10, 100, 30, 100, 1000, 1000)#0.95500000 0.01311106
#rayFisherVwInferenceExperiment(10, 100, 100, 100, 1000, 1000)#0.94700000 0.01416912
#rayFisherVwInferenceExperiment(30, 100, 30, 100, 1000, 1000)#0.94300000 0.01466301
#rayFisherVwInferenceExperiment(30, 100, 100, 100, 1000, 1000)#0.95900000 0.01254097
#rayFisherVwInferenceExperiment(100, 100, 100, 100, 1000, 1000)#0.95500000 0.01311106

# Experiments with equal means and sample sizes but differing kappas.
#rayFisherVwInferenceExperiment(3, 1, 3, 10, 1000, 2000)#0.91000000 0.01279844
#rayFisherVwInferenceExperiment(3, 1, 3, 100, 1000, 1000)#0.92700000 0.01645248
#rayFisherVwInferenceExperiment(3, 10, 3, 100, 1000, 1000)#0.94300000 0.01466301
#rayFisherVwInferenceExperiment(10, 1, 10, 10, 1000, 2000)#0.91500000 0.01247197
#rayFisherVwInferenceExperiment(10, 1, 10, 100, 1000, 2000)#0.90800000 0.01292563
#rayFisherVwInferenceExperiment(10, 10, 10, 100, 1000, 2000)#0.952000000 0.009559916
#rayFisherVwInferenceExperiment(30, 1, 30, 10, 1000, 2000)#0.9045000 0.0131438
#rayFisherVwInferenceExperiment(30, 1, 30, 100, 1000, 1000)#0.91700000 0.01744832
#rayFisherVwInferenceExperiment(30, 10, 30, 100, 1000, 1000)#0.96000000 0.01239355
#rayFisherVwInferenceExperiment(100, 1, 100, 10, 1000, 1000)#0.92300000 0.01686072
#rayFisherVwInferenceExperiment(100, 1, 100, 100, 1000, 1000)#0.91200000 0.01791714
#rayFisherVwInferenceExperiment(100, 10, 100, 100, 1000, 1000)#0.94500000 0.01441874

# Experiments on coverage rates using unequal means.
rayFisherVwInferenceUnequalExperiment <- function(n1, kappa1, n2, kappa2, numBoots, numTrials) {
  f <- function(i) {
    if (i %% 250 == 0)
      print(date())
    data1 <- rayFisher(rayUniform(), kappa1, n1)
    data2 <- rayFisher(rayUniform(), kappa2, n2)
    rayFisherVwInference(data1, data2, numBoots)
  }
  ps <- sapply(1:numTrials, f)
  pHat <- sum(ps > 0.05) / numTrials
  c(pHat, 2 * standardErrorProportion(numTrials, pHat))
}

# Experiments with unequal means, kappa1 = kappa2 = 1, varying sample size.
#rayFisherVwInferenceUnequalExperiment(3, 1, 3, 1, 1000, 1000)#0.89100000 0.01970979
#rayFisherVwInferenceUnequalExperiment(3, 1, 10, 1, 1000, 1000)#0.79500000 0.02553233
#rayFisherVwInferenceUnequalExperiment(3, 1, 30, 1, 1000, 1000)#0.73400000 0.02794595
#rayFisherVwInferenceUnequalExperiment(3, 1, 100, 1, 1000, 1000)#0.743000 0.027637
#rayFisherVwInferenceUnequalExperiment(10, 1, 10, 1, 1000, 1000)#0.57800000 0.03123562
#rayFisherVwInferenceUnequalExperiment(10, 1, 30, 1, 1000, 1000)#0.46700000 0.03155383
#rayFisherVwInferenceUnequalExperiment(10, 1, 100, 1, 1000, 1000)#0.40900000 0.03109463
#rayFisherVwInferenceUnequalExperiment(30, 1, 30, 1, 1000, 1000)#0.26000000 0.02774167
#rayFisherVwInferenceUnequalExperiment(30, 1, 100, 1, 1000, 1000)#0.16800000 0.02364538
#rayFisherVwInferenceUnequalExperiment(100, 1, 100, 1, 1000, 1000)#0.09000000 0.01809972

# Experiments with unequal means, kappa1 = kappa2 = 10, varying sample size.
#rayFisherVwInferenceUnequalExperiment(3, 10, 3, 10, 1000, 1000)#0.17900000 0.02424533
#rayFisherVwInferenceUnequalExperiment(3, 10, 10, 10, 1000, 1000)#0.10900000 0.01970979
#rayFisherVwInferenceUnequalExperiment(3, 10, 30, 10, 1000, 1000)#0.11400000 0.02010015
#rayFisherVwInferenceUnequalExperiment(3, 10, 100, 10, 1000, 1000)#0.1250000 0.0209165
#rayFisherVwInferenceUnequalExperiment(10, 10, 10, 10, 1000, 1000)#0.03500000 0.01162325
#rayFisherVwInferenceUnequalExperiment(10, 10, 30, 10, 1000, 1000)#0.023000000 0.009480717
#rayFisherVwInferenceUnequalExperiment(10, 10, 100, 10, 1000, 1000)#0.017000000 0.008175818
#rayFisherVwInferenceUnequalExperiment(30, 10, 30, 10, 1000, 1000)#0.013000000 0.007164077
#rayFisherVwInferenceUnequalExperiment(30, 10, 100, 10, 1000, 1000)#0.005000000 0.004460942
#rayFisherVwInferenceUnequalExperiment(100, 10, 100, 10, 1000, 1000)#0.003000000 0.003458902

# Experiments with unequal means, kappa1 = kappa2 = 100, varying sample size.
#rayFisherVwInferenceUnequalExperiment(3, 100, 3, 100, 1000, 1000)#0.015000000 0.007687652
#rayFisherVwInferenceUnequalExperiment(3, 100, 10, 100, 1000, 1000)#0.011000000 0.006596666
#rayFisherVwInferenceUnequalExperiment(3, 100, 30, 100, 1000, 1000)#0.012000000 0.006886509
#rayFisherVwInferenceUnequalExperiment(3, 100, 100, 100, 1000, 1000)#0.00700000 0.00527295
#rayFisherVwInferenceUnequalExperiment(10, 100, 10, 100, 1000, 1000)#0.00600000 0.00488426
#rayFisherVwInferenceUnequalExperiment(10, 100, 30, 100, 1000, 1000)#0.002000000 0.002825597
#rayFisherVwInferenceUnequalExperiment(10, 100, 100, 100, 1000, 1000)#0.002000000 0.002825597
#rayFisherVwInferenceUnequalExperiment(30, 100, 30, 100, 1000, 1000)#0.001000 0.001999
#rayFisherVwInferenceUnequalExperiment(30, 100, 100, 100, 1000, 1000)#0 0
#rayFisherVwInferenceUnequalExperiment(100, 100, 100, 100, 1000, 1000)#0 0

# Experiments with unequal means and kappas.
#rayFisherVwInferenceUnequalExperiment(3, 1, 3, 10, 1000, 1000)#0.73000000 0.02807846
#rayFisherVwInferenceUnequalExperiment(3, 1, 3, 100, 1000, 1000)#0.7540000 0.0272385
#rayFisherVwInferenceUnequalExperiment(3, 10, 3, 100, 1000, 1000)#0.11900000 0.02047818
#rayFisherVwInferenceUnequalExperiment(10, 1, 10, 10, 1000, 1000)#0.41900000 0.03120506
#rayFisherVwInferenceUnequalExperiment(10, 1, 10, 100, 1000, 1000)#0.41900000 0.03120506
#rayFisherVwInferenceUnequalExperiment(10, 10, 10, 100, 1000, 1000)#0.017000000 0.008175818
#rayFisherVwInferenceUnequalExperiment(30, 1, 30, 10, 1000, 1000)#0.13200000 0.02140804
#rayFisherVwInferenceUnequalExperiment(30, 1, 30, 100, 1000, 1000)#0.1350000 0.0216125
#rayFisherVwInferenceUnequalExperiment(30, 10, 30, 100, 1000, 1000)#0.002000000 0.002825597
#rayFisherVwInferenceUnequalExperiment(100, 1, 100, 10, 1000, 1000)#0.04000000 0.01239355
#rayFisherVwInferenceUnequalExperiment(100, 1, 100, 100, 1000, 1000)#0.046000 0.013249
#rayFisherVwInferenceUnequalExperiment(100, 10, 100, 100, 1000, 1000)#0.001000 0.001999



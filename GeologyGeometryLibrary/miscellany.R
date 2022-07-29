


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# These are miscellaneous functions to support the rest of the library. Some are preceded by detailed documentation. The others are not intended for use by end-users.



### DATA STRUCTURES ###

#' Higher-order function to apply n-nary function to n lists of arguments.
#' 
#' @param f An R function of n arguments.
#' @param ... n lists, where n is the same n as above. Denote them l1, ..., ln. They are assumed to have equal length l.
#' @return A list of length l, whose ith element is f(l1[[i]], ..., ln[[i]]).
thread <- function(f, ...)
  mapply(f, ..., SIMPLIFY=FALSE)

# Given a list of length n, returns a list of n-choose-2 lists of length 2.
allPairs <- function(set) {
  if (length(set) < 2)
    list()
  else
    c(lapply(set[2:length(set)], function(x) list(set[[1]], x)),
      allPairs(set[2:length(set)]))
}

# Given a list, returns a list consisting of sublists with each element of the list deleted.
listOmitting <- function(lst) {
  c(list(lst[2:length(lst)]),
    lapply(2:(length(lst) - 1), function(j) c(lst[1:(j - 1)], lst[(j + 1):length(lst)])),
    list(lst[1:(length(lst) - 1)]))
}

# Returns an empty data frame with the given number of rows.
# The implementation is stupid, but I can't find a clean answer for how to do this on the web.
emptyDataFrame <- function(nrow) {
  newFrame <- data.frame(dummy=replicate(nrow, NA))
  newFrame$dummy <- NULL
  newFrame
}

# Splits data frame into a list of rows, each with the correct names. I used to do this a lot, before I understood the data frame idiom. I still do it in ellipsoid tutorial 5dispersion.R.
listFromDataFrame <- function(df) {
  result <- list()
  for (i in 1:nrow(df)) {
    row <- list()
    for (name in names(df))
      row[[name]] <- df[[name]][[i]]
    result[[length(result) + 1]] <- row
  }
  result
}



### MISCELLANEOUS MATH ###

#' Conversion between degrees and radians.
#' 
#' Math works better in radians than in degrees. So R (like almost every programming language) operates in radians. So almost all of this library operates in radians. The exceptions are all labeled 'Deg', as in 'geoCartesianFromStrikeDipDeg', like a skull and crossbones on a bottle of poison. Anyway, to convert from degrees to radians, you multiply by this number. For example, 30 * degree == pi / 6. To convert from radians to degrees, you divide by this number. For example, 2 * pi / degree == 360. degree should really be defined in geology.R, but I don't want to call it geoDegree.
degree <- 2 * pi / 360

# epsilon is used in certain numerical methods, to establish stopping criteria.
epsilon <- .Machine$double.eps^0.25

#' General arithmetic mean operation for any vector space: scalars, vectors, matrices, etc.
#' 
#' @param xs A list. It is assumed that all items in the list are of the same type, and that that type supports the '+' operator.
#' @return An object of the same type as the elements in xs.
arithmeticMean <- function(xs) {
  Reduce("+", xs) / length(xs)
}

# Safe acos, for times when you know that you don't care about arguments slightly outside [-1, 1].
arcCos <- function(x) {
  acos(max(-1, min(1, x)))
}

arcSin <- function(x) {
  asin(max(-1, min(1, x)))
}

# Safe sqrt, for times when you know that you don't care about arguments slightly outside [0, infinity).
squareRoot <- function(x) {
  sqrt(max(0, x))
}

# Safe cube root, able to handle positive, zero, and negative arguments.
cubeRoot <- function(x) {
  if (x < 0)
    -(-x)^(1 / 3)
  else
    x^(1 / 3)
}

# Solutions of the quadratic polynomial a x^2 + b x + c == 0. Returns a numeric vector of 0, 1, 2, or 3 numbers. Typically there are 0 or 2 solutions. Occasionally the 2 solutions are identical. In sub-quadratic cases there could be 0, 1, or infinitely many solutions. The infinite case is signaled by returning 3 solutions.
realQuadraticSolutions <- function(aa, bb, cc) {
  if (aa == 0)
    if (bb == 0)
      if (cc == 0)
        # Infinitely many solutions; report three arbitrarily.
        c(0, 1, 2)
  else
    # No solutions.
    c()
  else
    # Typical linear case with one solution.
    c(-cc / bb)
  else {
    # Typical quadratic case with zero or two solutions.
    disc <- bb^2 - 4 * aa * cc
    if (disc < 0)
      c()
    else {
      root <- sqrt(disc)
      c((-bb + root) / (2 * aa), (-bb - root) / (2 * aa))
    }
  }
}

#' The standard error for a proportion.
#' 
#' As long as n pHat >= 10 and n (1 - pHat) >= 10, this is a reasonable approximation.
#' @param n A real number (positive integer).
#' @param pHat A real number (between 0 and 1).
#' @return A real number. The approximate 95% confidence region for the proportion is pHat plus or minus twice this number.
standardErrorProportion <- function(n, pHat) {
  sqrt(pHat * (1 - pHat) / n)
}

#' Trace of a matrix.
#' 
#' @param m A real nxn matrix.
#' @return A real number.
tr <- function(m) {
  sum(diag(m))
}

#' Cross product of three-dimensional Cartesian vectors.
#' 
#' @param v A 3D real vector.
#' @param w A 3D real vector.
#' @return A 3D real vector.
cross <- function(v, w) {
  c(v[[2]] * w[[3]] - v[[3]] * w[[2]], v[[3]] * w[[1]] - v[[1]] * w[[3]], v[[1]] * w[[2]] - v[[2]] * w[[1]])
}

#' Dot product of n-dimensional Cartesian vectors.
#' 
#' @param v A real vector.
#' @param w A real vector of the same dimension.
#' @return A real number.
dot <- function(v, w) {
  sum(v * w)
}

#' Antipodal map sending upper-hemisphere rays to the lower hemisphere.
#' 
#' In other words, negates any vector with positive [[3]]-component. Used in lower-hemisphere plots.
#' @param v A real 3D vector.
#' @return A real 3D vector.
lower <- function(v) {
  if (v[[3]] <= 0)
    v
  else
    -v
}

# Returns a new list of vectors, like vs with duplicates removed. Slow.
distinctVectors <- function(vs) {
  if (length(vs) <= 1)
    vs
  else {
    ws <- list(vs[[1]])
    for (i in 2:length(vs)) {
      # Test whether the ith vector matches any of those seen thus far.
      j <- 1
      while (j <= length(ws)) {
        diff <- vs[[i]] - ws[[j]]
        if (dot(diff, diff) < epsilon)
          # Don't include the ith vector.
          j <- length(ws) + 2
        else
          j <- j + 1
      }
      if (j == length(ws) + 1)
        # Include the ith vector.
        ws[[length(ws) + 1]] <- vs[[i]]
    }
    ws
  }
}

# The number interval [a, b] is mapped to the interval [0, 1]. Inputs outside [a, b] are clamped to [a, b].
scale <- function(y, a=0, b=1) {
  (y - a) / (b - a)
}

# The number interval range(ys) is mapped to the interval [0, 1].
scales <- function(ys) {
  sapply(ys, scale, min(ys), max(ys))
}

#' General machinery for permutation tests of regressions.
#' 
#' @param numPerms A real number (positive integer). The number of permutations to try.
#' @param f An R function. Must return a result list with fields $error, $minEigenvalue, $rSquared. Typically lineGeodesicRegression or a similar regression function.
#' @param xs A real vector. The values of the independent variable. Assumed to be the first argument passed to f.
#' @param ... Other arguments passed to f, after xs.
#' @return A real vector. The maximum length is numPerms. Often the length is less than numPerms, because the regression failed (as signaled by error != 0 or minEigenvalue <= 0). The user post-processes this vector to count how many of its entries exceed the original R^2 value. That fraction is a p-value for the null hypothesis that the regression result is meaningless.
permutedRSquareds <- function(numPerms, f, xs, ...) {
  print("fraction of permutations complete:")
  g <- function(i) {
    print(i / numPerms)
    newXs <- sample(xs, size=length(xs))
    regr <- f(xs=newXs, ...)
    c(regr$error, regr$minEigenvalue, regr$rSquared)
  }
  perms <- sapply(1:numPerms, g)
  perms[3,][perms[1,] == 0 & perms[2,] > 0]
}



### COORDINATE TRANSFORMATIONS ###

# Given (x, y, z), returns (rho, phi, theta), where 0 <= rho, 0 <= phi <= pi, and -pi <= theta <= pi.
sphericalFromCartesian <- function(xyz) {
  rho <- sqrt(dot(xyz, xyz))
  if (rho == 0)
    c(0, 0, 0)
  else {
    phi <- acos(xyz[[3]] / rho)
    rhoSinPhi <- rho * sin(phi)
    if (rhoSinPhi == 0)
      if (xyz[[3]] >= 0)
        c(xyz[[3]], 0, 0)
    else
      c(-xyz[[3]], pi, 0)
    else
      c(rho, phi, atan2(xyz[[2]], xyz[[1]]))
  }
}

# Given (rho, phi, theta), returns (x, y, z). Inverse to sphericalFromCartesian, up to periodicity.
cartesianFromSpherical <- function(rpt) {
  sinPhi <- sin(rpt[[2]])
  x <- rpt[[1]] * sinPhi * cos(rpt[[3]])
  y <- rpt[[1]] * sinPhi * sin(rpt[[3]])
  z <- rpt[[1]] * cos(rpt[[2]])
  c(x, y, z)
}

# Projection from a unit sphere out to the circumscribed unit cylinder parallel to the [[3]]-axis. Then drops the radial coordinate, because that is always 1. Area-preserving.
horizontalFromCartesian <- function(xyz) {
  if (xyz[[1]] == 0 && xyz[[2]] == 0)
    c(0, xyz[[3]])
  else
    c(atan2(xyz[[2]], xyz[[1]]), xyz[[3]])
}

# Inverse to horizontalFromCartesian.
cartesianFromHorizontal <- function(tz) {
  root <- sqrt(1 - tz[[2]]^2)
  c(cos(tz[[1]]) * root, sin(tz[[1]]) * root, tz[[2]])
}

# Assumes that xyz has length 1 and is not [0 0 1]^T.
equalAngleProjection <- function(xyz) {
  c(xyz[1] / (1 - xyz[[3]]), xyz[[2]] / (1 - xyz[[3]]))
}

# Assumes that xyz has length 1.
equalAreaProjection <- function(xyz) {
  normSquared <- xyz[[1]]^2 + xyz[[2]]^2
  if (normSquared == 0) {
    if (xyz[[2]] <= 0)
      # South pole.
      c(0, 0)
    else
      # North pole. This choice of point on the boundary is arbitrary.
      c(2, 0)
  }
  else {
    z <- min(1, max(-1, xyz[[3]]))
    sqrt(2 * (1 + z) / normSquared) * c(xyz[[1]], xyz[[2]])
  }
}



### ORDINARY DIFFERENTIAL EQUATIONS ###

#' Classical fourth-order Runge-Kutta method on a Euclidean space.
#' 
#' @param y0 A d-dimensional vector. The initial condition.
#' @param vel An R function that takes a time s and a d-dimensional vector y and returns a d-dimensional velocity vector v.
#' @param n A real number (positive integer). The number of steps to use.
#' @param correction A function to be applied to the vector after each step. Pass NULL to deactivate.
#' @return A d-dimensional vector. The final condition.
rungeKutta <- function(y0, vel, n, correction=NULL) {
  h <- 1 / n
  ys <- y0
  for (s in ((0:(n - 1)) / n)) {
    k1 <- vel(s, ys)
    y <- ys + k1 * h / 2
    if (!is.null(correction))
      y <- correction(y)
    k2 <- vel(s + h / 2, y)
    y <- ys + k2 * h / 2
    if (!is.null(correction))
      y <- correction(y)
    k3 <- vel(s + h / 2, y)
    y <- ys + k3 * h
    if (!is.null(correction))
      y <- correction(y)
    k4 <- vel(s + h, y)
    ys <- ys + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6
    if (!is.null(correction))
      ys <- correction(ys)
  }
  ys
}

#' Left-invariant fourth-order Runge-Kutta method on an arbitrary Lie group G.
#' 
#' From Munthe-Kaas (1998). Used to solve a differential equation of the form y-dot = y f, where y is a curve in G and f : G -> g is a Lie algebra-valued function describing the velocity of the curve when it is at y.
#' @param y0 An element of G. The initial condition at time t = 0.
#' @param f An R function that takes a Lie group element and returns a velocity vector in its tangent space.
#' @param n A real number (positive integer). The number of steps to use.
#' @param operation The group operation. The default suffices for matrix groups.
#' @param exponential The exponential map. The default suffices for matrix groups, but consider specialized cases such as rotExp.
#' @return An element of G. The final condition at time t = 1.
rungeKuttaLeft <- function(y0, f, n, operation=(function(g, h) {g %*% h}), exponential=expm) {
  h <- 1 / n
  ys <- y0
  for (s in ((0:(n - 1)) / n)) {
    k1 <- f(ys)
    k2 <- f(operation(ys, exponential(k1 * h / 2)))
    k3 <- f(operation(ys, exponential(k2 * h / 2 + (operation(k1, k2) - operation(k2, k1)) * h^2 / 24)))
    k4 <- f(operation(ys, exponential(k3 * h + (operation(k1, k3) - operation(k3, k1)) * h^2 / 6)))
    v <- (k1 + 2 * k2 + 2 * k3 + k4) * h / 6
    w <- (3 * k1 + 2 * k2 + 2 * k3 - k4) * h / 24
    ys <- operation(ys, exponential(v + (operation(w, v) - operation(v, w))))
  }
  ys
}

rungeKuttaRightDExpInv <- function(u, v, bracket) {
  v - bracket(u, v) / 2 + bracket(u, bracket(u, v)) / 24
}

#' Right-invariant fourth-order Runge-Kutta method on an arbitrary Lie group G.
#' 
#' From Munthe-Kaas (1999). Used to solve a differential equation of the form y-dot = f y, where y is a curve in G and f : R x G -> g is a Lie algebra-valued function describing the velocity of the curve when it is at time t and point y.
#' @param y0 An element of G. The initial condition at time t = 0.
#' @param f An R function that takes a time and a Lie group element and returns a velocity vector in its tangent space.
#' @param n A real number (positive integer). The number of steps to use.
#' @param operation The group operation. The default suffices for matrix groups.
#' @param exponential The exponential map. The default suffices for matrix groups, but consider specialized cases such as rotExp.
#' @param bracket The Lie bracket operation. The default suffices for matrix groups.
#' @return An element of G. The final condition at time t = 1.
rungeKuttaRight <- function(y0, f, n, operation=(function(g, h) {g %*% h}), exponential=expm, bracket=(function(u, v) {operation(u, v) - operation(v, u)})) {
  h <- 1 / n
  ys <- y0
  for (s in ((0:(n - 1)) / n)) {
    k1 <- f(h * s, ys)
    k2 <- rungeKuttaRightDExpInv(0.5 * h * k1, f(h * (s + 0.5), operation(exponential(0.5 * h * k1), ys)), bracket)
    k3 <- rungeKuttaRightDExpInv(0.5 * h * k2, f(h * (s + 0.5), operation(exponential(0.5 * h * k2), ys)), bracket)
    k4 <- rungeKuttaRightDExpInv(h * k3, f(h * (s + 1), operation(exponential(h * k3), ys)), bracket)
    ys <- operation(exponential((k1 + 2 * k2 + 2 * k3 + k4) * h / 6), ys)
  }
  ys
}

# Here is an example showing the Euclidean RK4 method to be a special case of the Lie group methods.
#y0 <- rnorm(3)
#l <- replicate(3, rnorm(3))
#l <- l - diag(c(1, 1, 1)) * tr(l) / 3
#rungeKutta(y0, function(s, y) {l %*% y}, 20, correction=NULL)
#rungeKuttaLeft(y0, function(y) {l %*% y}, 20, function(g, h) {g + h}, function(v) {v})
#rungeKuttaRight(y0, function(s, y) {l %*% y}, 20, function(g, h) {g + h}, function(v) {v})

# The first Lie group method can be used to solve d/dt(F) = L F, and the second Lie group method can be used to solve d/dt(F^-1) = F^-1 * -L, which should be equivalent. See the test after defLeftJeffery and defRightJeffery for an example.



### METRIC SPACES ###

#' The Euclidean distance between two d-dimensional points.
#' 
#' @param x A real vector.
#' @param y A real vector, of the same dimension as x.
#' @return A scalar. The distance between x and y.
euclideanDistance <- function(x, y) {
  sqrt(dot(x - y, x - y))
}

# Given a point x, a list of points seeds, and a distance function dist(x, seed, ...), returns the index of the seed to which x is closest.
voronoiCell <- function(x, seeds, dist, ...) {
  which.min(sapply(seeds, dist, x, ...))
}

#' Voronoi partition in an arbitrary metric space.
#' 
#' A clustering of n points in a metric space is a list of (usually non-empty) vectors of real numbers, such that all of the numbers are integers between 1 and n (inclusive) and each integer between 1 and n appears exactly once in the clustering. For example, list(c(5, 2), c(3, 1, 4)) is a clustering for n = 5.
#' @param xs A list of n points in the metric space.
#' @param seeds A list of k points in the metric space.
#' @param dist The metric on the metric space. An R function that takes two points as inputs and produces a non-negative real number as output. Additional arguments can be passed through ....
#' @return A list of k numeric vectors. The integers 1, 2, ..., n are partitioned among these lists. (Each of these integers appears once and only once among the lists.) The ith list is the indices of the xs that are in the Voronoi cell centered at the ith seed.
clusteringVoronoi <- function(xs, seeds, dist, ...) {
  clus <- replicate(length(seeds), c(), simplify=FALSE)
  for (i in 1:length(xs)) {
    j <- voronoiCell(xs[[i]], seeds, dist, ...)
    clus[[j]] <- c(clus[[j]], i)
  }
  clus
}

# DEPRECATED; USE clusteringVoronoi INSTEAD? Given n data points and m seed points, partitions the data points into m lists, which each list consists of the data points closest to one of the seed points. Some of the lists may be empty.
partitionVoronoi <- function(xs, seeds, dist, ...) {
  partition <- replicate(length(seeds), list(), simplify=FALSE)
  for (x in xs) {
    i <- which.min(sapply(seeds, dist, x, ...))
    partition[[i]][[length(partition[[i]]) + 1]] <- x
  }
  partition
}

#' Matrix of distances between points in a metric space.
#'
#' For an input list of length n, the algorithm uses O(n^2) calls to the distance function. So this can take a while.
#' @param xs A list of points in the metric space.
#' @param dist The metric on the metric space. An R function that takes two points as inputs and produces a non-negative real number as output. Additional arguments can be passed through ....
#' @param verbose A logical. Whether to print information about the progress of the computation.
#' @param ... Additional parameters to be passed to dist.
#' @return An n x n matrix, with only the top half filled in.
matrixOfDistances <- function(xs, dist, verbose=FALSE, ...) {
  # Compute the distances among the points, once and for all.
  n <- length(xs)
  d <- diag(replicate(n, 0))
  if (verbose)
    print("matrixOfDistances: verbose mode activated. Fraction complete...")
  for (i in 1:(n - 1)) {
    if (verbose)
      print(i / n)
    for (j in (i + 1):n)
      d[i, j] <- dist(xs[[i]], xs[[j]], ...)
  }
  # Copy the upper triangle into the lower triangle.
  #for (i in 2:n) {
  #  if (verbose)
  #    print(i / n)
  #  for (j in 1:(i - 1))
  #    d[i, j] <- d[j, i]
  #}
  d
}

#' DBSCAN clustering in an arbitrary metric space.
#'
#' A clustering of n points in a metric space is a list of (usually non-empty) vectors of real numbers, such that all of the numbers are integers between 1 and n (inclusive) and each integer between 1 and n appears exactly once in the clustering. For example, list(c(5, 2), c(3, 4, 1)) is a clustering for n = 5.
#' Intuitively, this algorithm forms clusters as contiguous patches of 'crowded' points. A point is 'crowded' if there are at least minPoints points within radius of it. You have to tune radius and minPoints to produce clusters that are meaningful to you.
#' @param dists An n x n matrix of pairwise distances between the points. Only the top half of the matrix is used. Typically obtained from matrixOfDistances.
#' @param radius A real number (positive). The radius about each point, that is considered that point's neighborhood.
#' @param minPoints A real number (positive integer). The minimum number of points needed, to consider a neighborhood crowded. Notice that the point itself is in its neighborhood. So no neighborhood ever has 0 neighbors.
#' @return A clustering. The final entry is a pseudo-cluster consisting of all of the points not included in any cluster. Unlike the real clusters, this pseudo-cluster may be empty.
clusteringDBSCAN <- function(dists, radius, minPoints) {
  # This is the only use of dists and radius.
  n <- ncol(dists)
  nbhd <- function(i) {
    Filter(function(j) {dists[min(i, j), max(i, j)] < radius}, 1:n)
  }
  # Prepare for the loop.
  clustering <- list()
  visited <- replicate(n, FALSE)
  clustered <- replicate(n, FALSE)
  i <- 1
  # Visit the points in order, although some will have been visited already when this outer loop reaches them.
  while (i <= n) {
    # Find the next unvisited point.
    while (i <= n && visited[[i]])
      i <- i + 1
    if (i <= n) {
      # i is the index of the first unvisited point. Visit it.
      visited[[i]] <- TRUE
      neighbors <- nbhd(i)
      if (length(neighbors) >= minPoints) {
        # Start a new cluster at point i.
        clus <- c(i)
        clustered[[i]] <- TRUE
        jj <- 1
        while (jj <= length(neighbors)) {
          j <- neighbors[[jj]]
          if (!visited[[j]]) {
            visited[[j]] <- TRUE
            newbors <- nbhd(j)
            if (length(newbors) >= minPoints)
              neighbors <- c(neighbors, newbors)
          }
          if (!clustered[[j]]) {
            clus[[length(clus) + 1]] <- j
            clustered[[j]] <- TRUE
          }
          jj <- jj + 1
        }
        # All points reachable from i are in clus.
        clustering[[length(clustering) + 1]] <- clus
      }
    }
    i <- i + 1
  }
  # Create a final pseudo-cluster from all of the not-yet-clustered noise points.
  noise <- Filter(function(i) {!clustered[[i]]}, 1:n)
  clustering[[length(clustering) + 1]] <- noise
  clustering
}

# Helper function for clusteringKMeans.
clusteringKMeansIntra <- function(clus, dd) {
  m <- outer(clus, clus, Vectorize(function(i, j) {dd[i, j]}))
  sum(m)
}

# Helper function for clusteringKMeans.
clusteringKMeansGamma <- function(dd, gamma) {
  # Compute the second term once for each cluster.
  f <- function(clus) {
    clusteringKMeansIntra(clus, dd) / length(clus)^2
  }
  seconds <- lapply(gamma, f)
  # Compute the best cluster l for each i.
  f <- function(l, i) {
    clus <- gamma[[l]]
    g <- function (r) {dd[i, clus[r]]}
    sum(sapply(1:length(clus), g)) * 2 / length(clus) - seconds[[l]]
  }
  new <- replicate(length(gamma), c(), simplify=FALSE)
  for (i in 1:(dim(dd)[1])) {
    l <- which.min(lapply(1:length(gamma), f, i))
    new[[l]] <- c(new[[l]], i)
  }
  new
}

#' k-means clustering in an arbitrary metric space, based on a seed clustering.
#'
#' A clustering of n points in a metric space is a list of (usually non-empty) vectors of real numbers, such that all of the numbers are integers between 1 and n (inclusive) and each integer between 1 and n appears exactly once in the clustering. For example, list(c(5, 2), c(3, 1, 4)) is a clustering for n = 5.
#' Intuitively, this algorithm forms clusters based on the idea that each point should be closer to its cluster's mean than to the mean of any other cluster. You have to supply the number k of clusters ahead of time. In fact, you have to supply a seed clustering. The algorithm improves that seed clustering until it stabilizes.
#' @param dists An n x n matrix of pairwise distances between the points. Only the top half of the matrix is used. Typically obtained from matrixOfDistances.
#' @param gamma0 A clustering. This seeds the algorithm. It can certainly affect the result. So you might want to run this function multiple times from different seeds.
#' @param numSteps A real number (non-negative integer). Bound on how many iterations to use.
#' @return A list consisting of gamma (a clustering), change (a non-negative real number), and numSteps (a non-negative integer). change measures how much the solution improved on the last iteration. numSteps is the number of iterations used.
clusteringKMeans <- function(dists, gamma0, numSteps=100) {
  # This isn't pretty, but let's make a full matrix of squared distances.
  dd <- dists^2
  for (i in 2:nrow(dists))
    for (j in 1:(i - 1))
      dd[i, j] <- dd[j, i]
  # Prepare for the loop.
  gamma <- gamma0
  wgd1 <- sum(sapply(gamma, clusteringKMeansIntra, dd))
  error <- epsilon + 1
  t <- 0
  # Update gamma until the change is small or the allotted time is exhausted.
  while (t < numSteps && error >= epsilon) {
    t <- t + 1
    wgd0 <- wgd1
    gamma <- clusteringKMeansGamma(dd, gamma)
    wgd1 <- sum(sapply(gamma, clusteringKMeansIntra, dd))
    error <- abs(wgd1 - wgd0)
  }
  list(gamma=gamma, change=error, numSteps=t)
}



### GRAPHICS ###

#' Assign a color to a number.
#' 
#' The number interval [a, b] is mapped to the color spectrum [red, magenta] (through blue and green). Inputs outside [a, b] are clamped to [a, b]. p is a distortion exponent. p = 1 was my original choice, but it under-uses the secondary colors. p = 2 makes the primaries a little under-used --- especially red. So I settled on the golden ratio.
#' @param y A real number. The number to be mapped to a color.
#' @param a A real number. The number that would be mapped to red.
#' @param b A real number. The number that would be mapped to magenta. Must not equal a.
#' @param p A real number (positive). The distortion exponent. I've never tried values outside 1 <= p <= 3.
#' @param opacity A real number (in [0, 1]). 0 is transparent, 1 is opaque. May not work in 3D plots.
#' @return Character. The color in '#rrggbb' format, suitable for passing to R graphics routines.
hue <- function(y, a=0, b=1, p=((1 + sqrt(5)) / 2), opacity=1) {
  # x is y rescaled and clamped to be between 0 and 1.
  x <- min(1, max(0, (y - a) / (b - a)))
  # z is x's fraction of the way from one benchmark to the next.
  z <- (x %% 0.2) / 0.2
  if (x < 0.2)
    # Interpolate from red to yellow by raising green.
    rgb(1, 1 - (1 - z)^p, 0, opacity)
  else if (x < 0.4)
    # Interpolate from yellow to green by lowering red.
    rgb(1 - z^p, 1, 0, opacity)
  else if (x < 0.6)
    # Interpolate from green to cyan by raising blue.
    rgb(0, 1, 1 - (1 - z)^p, opacity)
  else if (x < 0.8)
    # Interpolate from cyan to blue by lowering green.
    rgb(0, 1 - z^p, 1, opacity)
  else if (x < 1)
    # Interpolate from blue to magenta by raising red.
    rgb(1 - (1 - z)^p, 0, 1, opacity)
  else
    rgb(1, 0, 1, opacity)
}

#' Assign colors to a sequence of numbers.
#' 
#' The number interval range(ys) is mapped to the color spectrum [red, magenta] (through blue and green).
#' @param ys A vector of real numbers. The numbers to be mapped to colors.
#' @param p A real number (positive). The distortion exponent, as in the hue function.
#' @param opacity A real number (in [0, 1]). 0 is transparent, 1 is opaque. May not work in 3D plots.
#' @return A character vector, of the same length as ys. Each string is a color in '#rrggbb' format, suitable for passing to R graphics routines.
hues <- function(ys, p=((1 + sqrt(5)) / 2), opacity=1) {
  sapply(ys, hue, min(ys), max(ys), p, opacity)
}

#' Coloring the results of a cluster analysis.
#'
#' A clustering of n points in a metric space is a list of (usually non-empty) vectors of real numbers, such that all of the numbers are integers between 1 and n (inclusive) and each integer between 1 and n appears exactly once in the clustering. For example, list(c(5, 2), c(3, 1, 4)) is a clustering for n = 5. This function returns a vector of n colors, that can be used in plots of the cluster results.
#' @param clus A clustering of n points.
#' @param p A real number. The distortion exponent. See the hue function.
#' @return A character vector of length n. Exactly n colors are used, from red through green and blue to magenta.
clusteringHues <- function(clus, p=((1 + sqrt(5)) / 2)) {
  colors <- replicate(sum(sapply(clus, length)), "red")
  if (length(clus) >= 2)
    for (i in 2:length(clus))
      colors[clus[[i]]] <- hue(i, a=1, b=length(clus))
  colors 
}

#' Assign a grayscale shade to a number.
#' 
#' The given domain [a, b] is mapped to the given range [c, d], which should be a subset of [0, 1]. The range is then included into [0, 1], which is mapped into the color spectrum [black, white]. Inputs outside domain are clamped to domain. For example, if drawing on a white background, then maybe range should be c(0, 0.75).
#' @param y A real number. The number to be mapped to a shade.
#' @param domain A vector of two real numbers. The domain [a, b].
#' @param range A vector of two real numbers. The range [c, d]. Should obey 0 <= c < d <= 1.
#' @return Character. The color in '#rrggbb' format, suitable for passing to R graphics routines.
shade <- function(x, domain=c(0, 1), range=c(0, 1)) {
  y <- (x - domain[[1]]) * (range[[2]] - range[[1]]) / (domain[[2]] - domain[[1]]) + range[[1]]
  y <- min(range[[2]], max(range[[1]], y))
  rgb(y, y, y)
}

#' Assign grayscale shades to a sequence of numbers.
#' 
#' The number interval range(xs) is mapped to the grayscale range, where 0 is black and 1 is white.
#' @param xs A vector of real numbers. The numbers to be mapped to shades.
#' @param range A vector of two real numbers. The range [c, d]. Should obey 0 <= c < d <= 1. As in the shade function.
#' @return A character vector, of the same length as xs. Each string is a color in '#rrggbb' format, suitable for passing to R graphics routines.
shades <- function(xs, range=c(0, 1)) {
  sapply(xs, shade, domain=range(xs), range=range)
}

#' A pair of screenshots of an RGL window.
#' 
#' Once you have an RGL window with the image that you want, maximize it, to make it as large as possible on the screen. Keeping it open, return to R and invoke this function. Currently file names must end in 'png' and only PNG output is supported.
#' @param leftName Character. The name (or path) of the file, in which to store the first screenshot. NULL disables it.
#' @param rightName Character. The name (or path) of the file, in which to store the second screenshot. NULL disables it.
#' @return NULL
afterMaximizingWindow <- function(leftName=NULL, rightName=NULL, zoom=0.55) {
  if (!is.null(leftName)) {
    rgl.viewpoint(theta=45, phi=45)
    par3d(zoom=zoom)
    rgl.snapshot(leftName, fmt=substr(leftName, nchar(leftName) - 2, nchar(leftName)))
  }
  if (!is.null(rightName)) {
    rgl.viewpoint(theta=-45, phi=45)
    par3d(zoom=zoom)
    rgl.snapshot(rightName, fmt=substr(rightName, nchar(rightName) - 2, nchar(rightName)))
  }
}

# Shapes 15, 17, 19 are filled square, triangle, circle. Shapes 0, 2, 1 are unfilled.
# Traditionally lineations are squares, intermediates are triangles, and poles to foliation are circles.
# Styles could be 'solid', 'dotted', 'dashed', etc. See lty parameter in par.
plotEqualArea <- function(points=list(), curves=list(), colors=c("black"), shapes=c(19), styles=c("solid")) {
  # Set up the plot.
  r0 <- sqrt(2)
  r1 <- r0 + 0.1
  plot.new()
  plot.window(xlim=c(-r1, r1), ylim=c(-r1, r1),asp = 1)
  # Plot the points.
  if (length(points) >= 1) {
    xys <- sapply(points, function(v) equalAreaProjection(v))
    points(xys[1,], xys[2,], col=colors, pch=shapes)
  }
  # Plot the curves.
  if (length(curves) >= 1)
    for (i in 1:length(curves)) {
      xys <- sapply(curves[[i]], function(u) equalAreaProjection(u))
      lines(t(xys), lty=styles[[(i - 1) %% length(styles) + 1]])
    }
  # Plot the boundary circle.
  xys <- sapply(1:72, function(t) {
    theta <- t * 2.0 * pi / 72
    r0 * c(cos(theta), sin(theta))
  })
  polygon(xys[1,], xys[2,])
  # Plot some fringes, like on a hippie's leather jacket, man.
  lines(c(0.0, 0.0), c(r0, r1))
  lines(c(r0, r1), c(0.0, 0.0))
  lines(c(0.0, 0.0), -c(r0, r1))
  lines(-c(r0, r1), c(0.0, 0.0))
}

#' Infrastucture for 3D plots.
#'
#' Unless you are designing your own system of plotting, you probably do not want to call this function directly. Call rotEqualAnglePlot, etc. instead. See also rgl's plot3d.
#' @param radius A real number (non-negative). If not supplied, then it is deduced from the points, curves, and triangles.
#' @param points A list of 3D real vectors.
#' @param curves A list of lists of 3D real vectors.
#' @param triangles A list of length-3 lists of 3D real vectors.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param simplePoints A logical. Whether to plot points as points or as spheres. Spheres are better for conveying depth, but they may be slow for large data sets.
#' @param backgroundColor String (color). Color of background and fog.
#' @param curveColor A string (color). Color of curves.
#' @param curveWidth A real number (positive). Width of curves in pixels.
#' @param fogStyle A string, either 'none', 'linear', 'exp', or 'exp2'. The style of fog used to suggest depth. See help for rgl.bg.
#' @param pointSize A real number (positive). The size of the points or spheres used to depict points. If simplePoints, then measured in pixels with default 3. If not simplePoints, then measured in the same units as the radius of the plot, with default 0.02 * radius.
#' @return NULL.
plot3D <- function(radius=NULL, points=list(), curves=list(), triangles=list(), colors=c("white"), simplePoints=FALSE, backgroundColor="black", curveColor="white", curveWidth=1, fogStyle="linear", pointSize=NULL, axesColors=c("red", "green", "blue")) {
  # Initialize the window. The user will close it later.
  rgl.open()
  rgl.bg(color=backgroundColor, fogtype=fogStyle)
  # Figure out the radius and pointSize if necessary.
  if (is.null(radius)) {
    pointRadius <- 0
    curveRadius <- 0
    triangleRadius <- 0
    f <- function(pointList) {max(abs(simplify2array(pointList)))}
    if (length(points) >= 1)
      pointRadius <- f(points)
    if (length(curves) >= 1)
      curveRadius <- max(sapply(curves, f))
    if (length(triangles) >= 1)
      triangleRadius <- max(sapply(triangles, f))
    radius <- max(pointRadius, curveRadius, triangleRadius)
  }
  if (is.null(pointSize)) {
    if (simplePoints)
      pointSize <- 3
    else
      pointSize <- 0.02 * radius
  }
  # Draw the points.
  if (length(points) >= 1) {
    xs <- sapply(points, function(p) {p[1]})
    ys <- sapply(points, function(p) {p[2]})
    zs <- sapply(points, function(p) {p[3]})
    if (simplePoints)
      rgl.points(x=xs, y=ys, z=zs, color=colors, size=pointSize)
    else
      rgl.spheres(x=xs, y=ys, z=zs, radius=pointSize, color=colors, lit=FALSE)
  }
  # Draw the curves.
  for (cur in curves) {
    xs <- sapply(cur, function(p) {p[1]})
    ys <- sapply(cur, function(p) {p[2]})
    zs <- sapply(cur, function(p) {p[3]})
    rgl.linestrips(x=xs, y=ys, z=zs, color=curveColor, lwd=curveWidth)
  }
  # Draw the triangles.
  if (length(triangles) >= 1) {
    xs <- sapply(triangles, function(tri) {sapply(tri, function(p) {p[1]})})
    ys <- sapply(triangles, function(tri) {sapply(tri, function(p) {p[2]})})
    zs <- sapply(triangles, function(tri) {sapply(tri, function(p) {p[3]})})
    rgl.triangles(x=xs, y=ys, z=zs)#, alpha=1.0, back="cull")
  }
  # Draw the three coordinate axes.
  xs <- radius * c(0, 1, 1, 1, 1, 1, -1, -1, -1, -1)
  ys <- radius * c(0, 0, -0.1, 0.1, 0, 0, -0.1, 0.1, 0, 0)
  zs <- radius * c(0, 0, 0, 0, -0.1, 0.1, 0, 0, -0.1, 0.1)
  if (radius > 1) {
    xs <- c(xs, as.numeric(sapply(1:radius, function(x) c(x, x, x, x))))
    ys <- c(ys, replicate(floor(radius), c(-0.1, 0.1, 0, 0)))
    zs <- c(zs, replicate(floor(radius), c(0, 0, -0.1, 0.1)))
  }
  rgl.lines(x=xs, y=ys, z=zs, color=axesColors[[1]], lwd=3)
  rgl.lines(x=zs, y=xs, z=ys, color=axesColors[[2]], lwd=3)
  rgl.lines(x=ys, y=zs, z=xs, color=axesColors[[3]], lwd=3)
  NULL
}



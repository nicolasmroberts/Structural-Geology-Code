


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This file offers many functions for computing with rotations of 3D space, borrowed from all over the orientation statistics literature. In this file, 'rotation matrix' means 'real 3x3 matrix that is special orthogonal'. That's our most common representation of rotations, although that may change in a future release.



### CONVERSIONS AMONG REPRESENTATIONS OF ROTATIONS ###

# The radius of the equal-volume plot. Roughly 0.62.
rotEqualVolumeRadius <- (3 / (4 * pi))^(1 / 3)

#' Conversion from anti-symmetric matrix (infinitesimal rotation) to a vector of its non-redundant entries.
#'
#' @param w A rotation matrix.
#' @return A real 3D vector of length between 0 and pi, inclusive.
rotVectorFromAntisymmetric <- function(w) {
	c(w[3, 2], w[1, 3], w[2, 1])
}

#' Conversion from vector representation to anti-symmetric matrix (infinitesimal rotation).
#'
#' @param v A real 3D vector.
#' @return A rotation matrix.
rotAntisymmetricFromVector <- function(v) {
	matrix(c(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0), 3, 3)
}

#' Conversion from angle-axis representation to matrix.
#'
#' @param au A 4D real vector c(u1, u2, u3, a), with u unit length and a an angle in radians.
#' @return A rotation matrix.
rotMatrixFromAxisAngle <- function(ua) {
	m <- rotAntisymmetricFromVector(ua[1:3])
	diag(c(1, 1, 1)) + sin(ua[[4]]) * m + (1 - cos(ua[[4]])) * (m %*% m)
}

#' Conversion from matrix to angle-axis representation.
#'
#' @param r A rotation matrix.
#' @return A 4D real vector c(u1, u2, u3, a), with u unit length and a an angle in radians (between 0 and pi).
rotAxisAngleFromMatrix <- function(r) {
  cosine <- (tr(r) - 1) / 2
  if (cosine >= 1)
    # The angle is 0 and the axis doesn't matter.
    c(0, 0, -1, 0)
  else {
    # These are 2 * sin(a) >= 0 times the true values.
    u1 <- r[3, 2] - r[2, 3]
    u2 <- r[1, 3] - r[3, 1]
    u3 <- r[2, 1] - r[1, 2]
    normU <- squareRoot(u1^2 + u2^2 + u3^2)
    if (normU != 0)
      c(u1 / normU, u2 / normU, u3 / normU, arcCos(cosine))
    else if (cosine > 0)
      # The angle is 0 and the axis doesn't matter.
      c(0, 0, -1, 0)
    else {
      # The angle is pi and the axis is unstable. Guess the absolute values of u1, u2, u3.
      u1 <- squareRoot((r[1, 1] + 1) / 2)
      u2 <- squareRoot((r[2, 2] + 1) / 2)
      u3 <- squareRoot((r[3, 3] + 1) / 2)
      # We can choose signs based on Rij = 2 ui uj when i != j.
      if (u1 != 0) {
        # Pick u1 > 0.
        if (r[1, 2] < 0)
          u2 <- -u2
        if (r[1, 3] < 0)
          u3 <- -u3
      } else if (u2 != 0) {
        # Pick u2 > 0.
        if (r[2, 1] < 0)
          u1 <- -u1
        if (r[2, 3] < 0)
          u3 <- -u3
      } else {
        # Pick u3 > 0.
        if (r[3, 1] < 0)
          u1 <- -u1
        if (r[3, 2] < 0)
          u2 <- -u2
      }
      c(u1, u2, u3, pi)
    }
  }
}

rotAxisAngleFromMatrixDeprecated <- function(r) {
	cosine <- (tr(r) - 1) / 2
	if (cosine >= 1)
		c(0, 0, -1, 0)
	else {
		u <- squareRoot(1 + (r[1, 1] - 1) / (1 - cosine))
		if (r[3, 2] < r[2, 3])
			u <- -u
		v <- squareRoot(1 + (r[2, 2] - 1) / (1 - cosine))
		if (r[1, 3] < r[3, 1])
			v <- -v
		w <- squareRoot(1 + (r[3, 3] - 1) / (1 - cosine))
		if (r[2, 1] < r[1, 2])
			w <- -w
		nrm <- sqrt(u^2 + v^2 + w^2)
		a <- arcCos(cosine)
		c(u / nrm, v / nrm, w / nrm, a)
	}
}

#' Matrix exponentiation of infinitesimal rotation to finite rotation.
#'
#' @param w A 3x3 real matrix (anti-symmetric).
#' @return A rotation matrix.
rotExp <- function(w) {
	th <- sqrt(tr(crossprod(w, w)) / 2)
	if (th == 0)
    diag(c(1, 1, 1))
	else
    diag(c(1, 1, 1)) + (sin(th) / th) * w + ((1 - cos(th)) / th^2) * (w %*% w)
}

#' Matrix logarithm to produce infinitesimal from finite rotation.
#'
#' @param r A 3x3 rotation matrix.
#' @return A 3x3 real matrix (anti-symmetric). The principal logarithm.
rotLog <- function(r) {
	ua <- rotAxisAngleFromMatrix(r)
	ua[[4]] * rotAntisymmetricFromVector(ua[1:3])
}

#' Maps a tangent space into the space of rotations.
#' 
#' Converts tangent vector v into rotation matrix R = exp(V) C.
#' @param v A 3D real vector.
#' @param center A rotation matrix.
#' @return A rotation matrix.
rotMatrixFromRightTangent <- function(v, center) {
  nrm <- sqrt(dot(v, v))
  if (nrm == 0)
    center
  else
    rotMatrixFromAngleAxis(c(v / nrm, nrm)) %*% center
}

#' Maps the space of rotations into one of its tangent spaces.
#' 
#' Converts rotation matrix R = exp(V) C into vector v. Appropriate only if the two given rotations are close to each other.
#' @param r A rotation matrix.
#' @param center A rotation matrix.
#' @return A 3D real vector.
rotRightTangentFromMatrix <- function(r, center) {
  ua <- rotAxisAngleFromMatrix(r %*% t(center))
  ua[1:3] * ua[[4]]
}

#' Maps a tangent space into the space of rotations.
#' 
#' Converts tangent vector v into rotation matrix R = C exp(V).
#' @param v A 3D real vector.
#' @param center A rotation matrix.
#' @return A 3x3 rotation matrix.
rotMatrixFromLeftTangent <- function(v, center) {
  nrm <- sqrt(dot(v, v))
  if (nrm == 0)
    center
  else
    center %*% rotMatrixFromAxisAngle(c(v / nrm, nrm))
}

#' Maps the space of rotations into one of its tangent spaces.
#' 
#' Converts rotation matrix R = exp(V) C into vector v. Appropriate only if the two given rotations are close to each other.
#' @param r A rotation matrix.
#' @param center A rotation matrix.
#' @return A 3D real vector.
rotLeftTangentFromMatrix <- function(r, center) {
  ua <- rotAxisAngleFromMatrix(t(center) %*% r)
  ua[1:3] * ua[[4]]
}

#' Rotation matrix about the x-axis.
#'
#' @param a A real number (angle in radians).
#' @return A rotation matrix.
rotMatrixAboutX <- function(a) {
	cosine <- cos(a)
	sine <- sin(a)
	matrix(c(1, 0, 0, 0, cosine, sine, 0, -sine, cosine), 3, 3)
}

#' Rotation matrix about the y-axis.
#'
#' @param a A real number (angle in radians).
#' @return A rotation matrix.
rotMatrixAboutY <- function(a) {
	cosine <- cos(a)
	sine <- sin(a)
	matrix(c(cosine, 0, -sine, 0, 1, 0, sine, 0, cosine), 3, 3)
}

#' Rotation matrix about the z-axis.
#'
#' @param a A real number (angle in radians).
#' @return A rotation matrix.
rotMatrixAboutZ <- function(a) {
	cosine <- cos(a)
	sine <- sin(a)
	matrix(c(cosine, sine, 0, -sine, cosine, 0, 0, 0, 1), 3, 3)
}

#' Rotation matrix from xzx Euler angles.
#'
#' The resulting matrix represents rotation about the global x-axis through the angle cba[[3]], followed by rotation about the global z-axis through cba[[2]], followed by rotation about the global x-axis through cba[[1]].
#' @param cba A 3D real vector. The Euler angles in radians.
#' @return A rotation matrix.
rotMatrixFromXZXAngles <- function(cba) {
  rotMatrixAboutX(cba[[1]]) %*% rotMatrixAboutZ(cba[[2]]) %*% rotMatrixAboutX(cba[[3]])
}

#' Extraction of xzx Euler angles from rotation matrix.
#'
#' A vector cba is produced, such that the matrix represents rotation about the global x-axis through the angle cba[[3]], followed by rotation about the global z-axis through cba[[2]], followed by rotation about the global x-axis through cba[[1]].
#' @param r A rotation matrix.
#' @return A 3D real vector. The Euler angles in radians. The middle angle is always in [0, pi].
rotXZXAnglesFromMatrix <- function(r) {
  bb <- arcCos(r[1, 1])
  aa <- atan2(r[1, 3], -r[1, 2])
  cc <- atan2(r[3, 1], r[2, 1])
  c(cc, bb, aa)
}

#' Conversion of angle-axis representation to equal-volume representation.
#'
#' @param ua A 4D real vector c(u1, u2, u3, a), with u unit length and a an angle in radians.
#' @return A 3D real vector.
rotEqualVolumeFromAxisAngle <- function(ua) {
	rhoCubed <- (ua[[4]] - sin(ua[[4]])) * 3 / (4 * pi^2)
	rhoCubed^(1 / 3) * ua[1:3]
}

#' Solution of x - sin(x) == d, within tolerance of epsilon.
#'
#' Helper function for rotAxisAngleFromEqualVolume, etc.
#' @param d A real number.
#' @return A real number, hopefully in [0, pi].
rotXMinusSinXSolution <- function(d) {
	# On [0, pi], x^3 / 9 is close to x - sin(x). Use that to choose the seed.
	x0 <- (9.0 * d)^(1.0 / 3.0)
	# Proceed by custom Newton's method. (See also optim, nlm, nlminb.)
	x1 <- x0 - (x0 - sin(x0) - d) / (1.0 - cos(x0))
	while (abs(x1 - x0) > epsilon) {
		x0 <- x1
		x1 <- x0 - (x0 - sin(x0) - d) / (1.0 - cos(x0))
	}
	x1
}

#' Conversion of equal-volume representation to angle-axis representation.
#'
#' @param v A 3D real vector.
#' @return A 4D real vector c(u1, u2, u3, a), with u unit length and a an angle in radians.
rotAxisAngleFromEqualVolume <- function(v) {
	rho <- sqrt(dot(v, v))
	if (rho == 0)
    c(0, 0, -1, 0)
	else {
		a <- rotXMinusSinXSolution(rho^3 * 4 * pi^2 / 3)
    c(v / rho, a)
	}
}

#' Conversion of rotation matrix to its equal-volume representation.
#'
#' @param r A rotation matrix.
#' @return A 3D real vector.
rotEqualVolumeFromMatrix <- function(r) {
	rotEqualVolumeFromAxisAngle(rotAxisAngleFromMatrix(r))
}

#' Conversion of equal-volume representation to rotation matrix.
#'
#' @param v A 3D real vector.
#' @return A rotation matrix.
rotMatrixFromEqualVolume <- function(v) {
	rotMatrixFromAxisAngle(rotAxisAngleFromEqualVolume(v))
}

#' Conversion of rotation matrix to its equal-angle representation.
#'
#' @param r A rotation matrix.
#' @return A 3D real vector.
rotEqualAngleFromMatrix <- function(r) {
  ua <- rotAxisAngleFromMatrix(r)
  ua[1:3] * tan(ua[[4]] / 4)
}

#' Conversion of equal-angle representation to rotation matrix.
#'
#' @param v A 3D real vector.
#' @return A rotation matrix.
rotMatrixFromEqualAngle <- function(v) {
  nrm <- sqrt(dot(v, v))
  if (nrm == 0)
    diag(c(1, 1, 1))
  else
    rotMatrixFromAxisAngle(c(v / nrm, 4 * atan(nrm)))
}

#' Conversion of rotation matrix to its Rodrigues representation.
#'
#' @param r A rotation matrix.
#' @return A 3D real vector.
rotRodriguesFromMatrix <- function(r) {
  ua <- rotAxisAngleFromMatrix(r)
  ua[1:3] * tan(ua[[4]] / 2)
}

#' Conversion of Rodrigues representation to rotation matrix.
#'
#' @param v A 3D real vector.
#' @return A rotation matrix.
rotMatrixFromRodrigues <- function(v) {
  nrm <- sqrt(dot(v, v))
  if (nrm == 0)
    diag(c(1, 1, 1))
  else
    rotMatrixFromAxisAngle(c(v / nrm, 2 * atan(nrm)))
}

#' Conversion of axis-angle representation to quaternion.
#'
#' @param ua A 4D real vector. Unit vector u followed by angle a in radians.
#' @return A 4D real vector (unit length).
rotQuaternionFromAxisAngle <- function(ua) {
  c(cos(ua[[4]] / 2), sin(ua[[4]] / 2) * ua[1:3])
}

#' Conversion of quaternion representation to axis-angle.
#'
#' @param q A 4D real vector (unit length).
#' @return A 4D real vector. Unit vector u followed by angle a in radians.
rotAxisAngleFromQuaternion <- function(q) {
	a <- 2 * arcCos(q[[1]])
	sine <- sin(a / 2)
	v <- q[2:4]
	if (sine == 0)
		if (abs(v[[1]]) > abs(v[[2]]) && abs(v[[1]]) > abs(v[[3]]))
			if (v[[1]] > 0)
        c(1, 0, 0, a)
			else
        c(-1, 0, 0, a)
		else if (abs(v[[2]]) > abs(v[[1]]) && abs(v[[2]]) > abs(v[[3]]))
			if (v[[2]] > 0)
        c(0, 1, 0, a)
			else
        c(0, -1, 0, a)
		else
			if (v[[3]] > 0)
        c(0, 0, 1, a)
			else
        c(0, 0, -1, a)
	else
    c(v / sine, a)
}

#' Conversion of quaternion representation to its rotation matrix.
#'
#' @param q A 4D real vector (unit length).
#' @return A rotation matrix.
rotMatrixFromQuaternion <- function(q) {
	r11 <- q[[1]]^2 + q[[2]]^2 - q[[3]]^2 - q[[4]]^2
	r12 <- 2 * (q[[2]] * q[[3]] - q[[1]] * q[[4]])
	r13 <- 2 * (q[[1]] * q[[3]] + q[[2]] * q[[4]])
	r21 <- 2 * (q[[1]] * q[[4]] + q[[2]] * q[[3]])
	r22 <- q[[1]]^2 + q[[3]]^2 - q[[2]]^2 - q[[4]]^2
	r23 <- 2 * (q[[3]] * q[[4]] - q[[1]] * q[[2]])
	r31 <- 2 * (q[[2]] * q[[4]] - q[[1]] * q[[3]])
	r32 <- 2 * (q[[1]] * q[[2]] + q[[3]] * q[[4]])
	r33 <- q[[1]]^2 + q[[4]]^2 - q[[2]]^2 - q[[3]]^2
	matrix(c(r11, r21, r31, r12, r22, r32, r13, r23, r33), 3, 3)
}

#' Conversion of rotation matrix to its quaternion representation.
#'
#' @param r A rotation matrix.
#' @return A 4D real vector (unit length).
rotQuaternionFromMatrix <- function(r) {
	rotQuaternionFromAxisAngle(rotAxisAngleFromMatrix(r))
}



### GEODESIC METHODS ###

#' The distance between two rotations as points in SO(3).
#'
#' @param r A rotation matrix.
#' @param q A rotation matrix.
#' @return A real number (in the interval [0, pi]).
rotDistance <- function(r, q) {
  arcCos((tr(crossprod(r, q)) - 1) / 2)
}

#' Diameter of a set of rotations.
#' 
#' @param rs A list of rotation matrices.
#' @return A real number (in the interval [0, pi]).
rotDiameter <- function(rs) {
  f <- function(i) {
    max(sapply(1:(i - 1), function(j) rotDistance(rs[[i]], rs[[j]])))
  }
  max(sapply(2:(length(rs)), f))
}

#' The smallest rotation matrix R such that R u = v.
#' 
#' @param u A ray (unit real 3D vector).
#' @param v A ray.
#' @return A rotation matrix.
rotSmallestRotationFromTwoRays <- function(u, v) {
  axis <- rayNormalized(cross(u, v))
  angle <- arcCos(dot(u, v))
  rotMatrixFromAxisAngle(c(axis, angle))
}

#' The smallest rotation matrix R such that R u = v or R u = -v.
#' 
#' @param u A line (unit real 3D vector).
#' @param v A line.
#' @return A rotation matrix.
rotSmallestRotationFromTwoLines <- function(u, v) {
  if (dot(u, v) < 0)
    rotSmallestRotationFromTwoRays(u, -v)
  else
    rotSmallestRotationFromTwoRays(u, v)
}

#' The minimum distance between two elements of a set of rotations.
#'
#' @param rs A list of rotation matrices.
#' @return A real number (in the interval [0, pi]).
rotSeparation <- function(rs) {
  f <- function(i) {
    min(sapply(1:(i - 1), function(j) rotDistance(rs[[i]], rs[[j]])))
  }
  min(sapply(2:(length(rs)), f))
}

#' The Frechet (geodesic L^2) variance of a set of rotations about a given rotation.
#'
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix.
#' @return A real number (between 0 and pi^2 / 2).
rotVariance <- function(rs, center) {
  dists <- sapply(rs, rotDistance, center)
  sum(dists^2) / (2 * length(rs))
}

#' The Frechet (geodesic L^2) mean of a set of rotations.
#'
#' An interative algorithm for computing the Frechet mean --- the rotation that minimizes the Frechet variance. The iterations continue until error squared of epsilon is achieved or numSteps iterations have been used. Try multiple seeds, to improve your chances of finding the global optimum.
#' @param rs A list of rotation matrices.
#' @param numSeeds A real number (positive integer). How many rs to try as seeds.
#' @param numSteps A real number (positive integer). Bound on how many iterations to use.
#' @return A list consisting of $mean (a special orthogonal real 3x3 matrix), $variance (a real number), $changeSquared (a real number), and $numSteps (a non-negative integer). changeSquared is the square of the size of the final step. numSteps is the number of iterations used.
rotMeanVariance <- function(rs, numSeeds=1, numSteps=100) {
  seeds <- sample(rs, numSeeds)
  # No variance is ever as large as 5.
  best <- c(5)
  for (seed in seeds) {
    rBar <- seed
    changeSquared <- epsilon + 1.0
    k <- 0
    while (changeSquared >= epsilon && k < numSteps) {
      w <- arithmeticMean(lapply(rs, function(r) rotLog(crossprod(rBar, r))))
      rBar <- rBar %*% rotExp(w)
      changeSquared <- tr(crossprod(w, w))
      k <- k + 1
    }
    var <- rotVariance(rs, rBar)
    if (var < best[[1]])
      best <- list(var, rBar, changeSquared, k)
  }
  list(variance=best[[1]], mean=best[[2]], changeSquared=best[[3]], numSteps=best[[4]])
}

#' The Frechet (geodesic L^2) mean. Convenience shortcut for rotMeanVariance.
#' 
#' @param rs A list of rotation matrices.
#' @param numSeeds A real number (positive integer). How many rs to try as seeds.
#' @param numSteps A real number (positive integer). Bound on how many iterations to use.
#' @return A rotation matrix.
rotFrechetMean <- function(rs, numSeeds=1, numSteps=100) {
  rotMeanVariance(rs, numSeeds=numSeeds, numSteps=numSteps)$mean
}

#' Point on a geodesic that is closest to a given point.
#' 
#' Returns the point on the geodesic (exp a M) B that is closest to Q.
#' @param q A rotation matrix.
#' @param m A 3x3 real matrix (antisymmetric).
#' @param b A rotation matrix.
#' @return A rotation matrix.
rotNearestPointOnRightGeodesic <- function(q, m, b) {
  v <- rotVectorFromAntisymmetric(m)
  vNorm <- sqrt(dot(v, v))
  if (vNorm == 0)
    b
  else {
    u <- m / vNorm
    uBQT <- u %*% b %*% t(q)
    alpha <- atan2(-tr(uBQT), tr(u %*% uBQT))
    r <- rotExp(alpha * u) %*% b
    s <- rotExp((alpha + pi) * u) %*% b
    if (rotDistance(q, r) < rotDistance(q, s))
      r
    else
      s
  }
}

#' Point on a geodesic that is closest to a given point.
#' 
#' Returns the point on the geodesic B (exp a M) that is closest to Q.
#' @param q A rotation matrix.
#' @param m A 3x3 real matrix (antisymmetric).
#' @param b A rotation matrix.
#' @return A rotation matrix.
rotNearestPointOnLeftGeodesic <- function(q, m, b) {
  v <- rotVectorFromAntisymmetric(m)
  vNorm <- sqrt(dot(v, v))
  if (vNorm == 0)
    b
  else {
    u <- m / vNorm
    qTBU <- t(q) %*% b %*% u
    alpha <- atan2(-tr(qTBU), tr(qTBU %*% u))
    r <- b %*% rotExp(alpha * u)
    s <- b %*% rotExp((alpha + pi) * u)
    if (rotDistance(q, r) < rotDistance(q, s))
      r
    else
      s
  }
}
# Here's some test code.
#b <- rotUniform()
#u <- rayUniform()
#rs <- lapply(0:360, function(i) b %*% rotMatrixFromAxisAngle(c(u, i * degree)))
#m <- rotAntisymmetricFromVector(1.34 * u) # arbitrary number
#q <- rotUniform()
#s <- rotNearestPointOnLeftGeodesic(q, m, b)
#min(sapply(rs, function(r) rotDistance(q, r)))
#rotDistance(q, s)
#rotEqualVolumePlot(list(b, s, q), list(rs, rotGeodesicPoints(q, s, 100)))

#' Points evenly spaced on a geodesic from one rotation to another.
#'
#' Doesn't work well if rotations are pi away from each other. Give it an intermediate point, to help it out.
#' @param r A rotation matrix.
#' @param q A rotation matrix.
#' @param numSteps A real number (positive integer).
#' @return A list of rotation matrices, of length numSteps + 1. The first one is r and the last one is q.
rotGeodesicPoints <- function(r, q, numSteps=10) {
  ua <- rotAxisAngleFromMatrix(r %*% t(q))
  lapply(0:numSteps, function(i) rotMatrixFromAxisAngle(c(ua[1:3], ua[[4]] * i / numSteps)) %*% q)
}



### PROJECTED ARITHMETIC MEAN AND RELATED COMPUTATIONS ###

#' Projects matrices near SO(3) onto SO(3).
#'
#' @param m A real 3x3 matrix. Presumed to be nearly special orthogonal.
#' @return A rotation matrix.
rotProjectedMatrix <- function(m) {
  # Q D^(1 / 2) Q^T = (M^T M)^(1 / 2).
	valsvecs <- eigen(crossprod(m, m), symmetric=TRUE)
	q <- valsvecs$vectors
	dSqrtInv <- valsvecs$values^(-1 / 2)
	# projection = M Q D^(-1 / 2) Q^T.
	m %*% q %*% diag(dSqrtInv) %*% t(q)
}

#' The projected arithmetic mean of a set of rotations.
#'
#' This function is appropriate only if the rotations are already known to be clustered about a central tendency (rather than girdled, say). In this case the projected arithmetic mean equals the MLE of the matrix Fisher mean (rotFisherMLE) and the quaternionic line mean (rotMeanScatter). If the data are not clustered, then the projected arithmetic mean may equal the negation of the quaternionic line mean (so orthogonal but not special orthogonal).
#' @param rs A list of rotation matrices.
#' @return A rotation matrix.
rotProjectedMean <- function(rs) {
  rotProjectedMatrix(arithmeticMean(rs))
}

#' Mean and dispersion of a sample, computed via line treatment of quaternions.
#' 
#' See Prentice (1986), p. 218.
#' @param rs A list of rotation matrices.
#' @return A list consisting of $values (real 4-vector, non-negative, descending order, sum to 1) and $rotations (list of four special orthogonal real 3x3 matrices). If val1 + val4 > 0.5, then the sample is bipolar and $rotations[[1]] equals the projected arithmetic mean and the MLE of the matrix Fisher mean. If val1 + val4 < 0.5, then the sample is equatorial and $rotations[[4]] equals the MLE of the matrix Fisher mean.
rotMeanScatter <- function(rs) {
  qs <- lapply(rs, rotQuaternionFromMatrix)
  ts <- lapply(qs, function(q) {outer(q, q)})
  tMatrix <- arithmeticMean(ts)
  eig <- eigen(tMatrix, symmetric=TRUE)
  rots <- lapply(1:4, function(j) rotMatrixFromQuaternion(eig$vectors[,j]))
  list(values=eig$values, rotations=rots)
}

# Given four rays, returns the rotation matrix R that takes a0 to a1 and, as well as possible, also b0 to b1. 
rotMatrixFromFourRaysAsymmetric <- function(a0, b0, a1, b1) {
  c0 <- rayOrthogonalProjection(a0, b0)
  c1 <- rayOrthogonalProjection(a1, b1)
  first <- cbind(a1, c1, cross(a1, c1))
  second <- cbind(a0, c0, cross(a0, c0))
  rotProjectedMatrix(first %*% t(second))
}

# Test code. R takes a0 to a1. Typically it does not take b0 to b1. But no other rotation S taking a0 to a1 produces a smaller angle between S b0 and b1.
rotMatrixFromFourRaysAsymmetricTest <- function() {
  a0 <- rayUniform()
  a1 <- rayUniform()
  b0 <- rayUniform()
  b1 <- rayUniform()
  r <- rotMatrixFromFourRaysAsymmetric(a0, b0, a1, b1)
  #print(det(r))
  #print(r %*% t(r))
  print(as.numeric(r %*% a0))
  print(a1)
  print(as.numeric(r %*% b0))
  print(b1)
  goal <- dot(as.numeric(r %*% b0), b1)
  print(goal)
  print(range(sapply(
    seq(from=0, to=(2 * pi), by=0.1), 
    function(s) dot(as.numeric(rotMatrixFromAxisAngle(c(a1, s)) %*% r %*% b0), b1))))
}

# Returns the rotation matrix R that approximately takes a0 to a1 and b0 to b1. If the angle between a0 and b0 does not match the angle between a1 and b1, then symmetrically error-corrects it (meaning, returns the rotation that takes a0 to a2 and b0 to b2, where the angle between a1 and a2 equals the angle between b1 and b2).
rotMatrixFromFourRays <- function(a0, b0, a1, b1) {
  c0 <- rayNormalized(a0 + b0)
  c1 <- rayNormalized(a1 + b1)
  rotMatrixFromFourRaysAsymmetric(c0, b0, c1, b1)
}

# Test code.
rotMatrixFromFourRaysTest <- function() {
  a0 <- rayUniform()
  a1 <- rayUniform()
  b0 <- rayUniform()
  b1 <- rayUniform()
  r <- rotMatrixFromFourRays(a0, b0, a1, b1)
  print(dot(as.numeric(r %*% a0), a1))
  print(dot(as.numeric(r %*% b0), b1))
}

# a0 and a1 are rays. b0 and b1 are lines. Returns the rotation that takes a0 to a1 and b0 to b1 approximately, in the same manner as rotMatrixFromFourRays. Actually there are two such rotations; returns the smaller.
rotMatrixFromTwoRaysTwoLines <- function(a0, b0, a1, b1) {
  r <- rotMatrixFromFourRays(a0, b0, a1, b1)
  s <- rotMatrixFromFourRays(a0, b0, a1, -b1)
  # dist(R, I) < dist(S, I) iff tr R > tr S.
  if (tr(r) > tr(s))
    r
  else
    s
}

# Helper function for rotationSVD. Given permutation P of 1:n,
# returns matrix M of 0s and 1s such that for all V, (M V)[i] == V[P[i]].
rotOrthogonalFromPermutation <- function(p) {
	n <- length(p)
	m <- matrix(0, n, n)
	for (i in p)
		m[i, p[i]] <- 1
	m
}

#' Signed singular value decomposition, using only rotations but allowing one negative singular value.
#'
#' @param m A real 3x3 matrix.
#' @return A list of three 3x3 real matrices $u, $d, $v, such that M = U D V^T. U and V are rotation matrices and D is diagonal with |D_11| >= D_22 >= D_33 >= 0. If M is an arithmetic mean of rotation matrices, then furthermore 1 >= |D_11|.
rotSingularValueDecomposition <- function(m) {
	# M = U D V^T.
	duv <- svd(m)
	u <- duv$u
	vT <- t(duv$v)
	# Ensure that the singular values are in decreasing order.
	p <- rotOrthogonalFromPermutation(order(duv$d, decreasing=TRUE))
	u <- u %*% t(p)
	d <- p %*% diag(duv$d) %*% t(p)
	vT <- p %*% vT
	# Ensure that det U = det V = 1.
	if (det(u) < 0.0) {
		u[,1] <- -u[,1]
		d[1, 1] <- -d[1, 1]
	}
	if (det(vT) < 0.0) {
		vT[1,] <- -vT[1,]
		d[1, 1] <- -d[1, 1]
	}
	list(u=u, d=d, v=t(vT))
}

#' P-value function for whether the data come from a distribution with the hypothesized mean.
#'
#' Based on sampling theory of the moment of inertia matrix (Prentice, 1986).
#' @param rs A list of rotation matrices.
#' @return An R function from {rotation matrices} to {real numbers union NULL}. If the data are unsuitable (because they are not bipolar), then NULL is always output. Otherwise, for any given hypothesized mean, this function produces the p-value for that hypothesis.
rotPrenticeInference <- function(rs) {
  # Convert to quaternions in unusual order.
  qs <- lapply(rs, rotQuaternionFromMatrix)
  qs <- lapply(qs, function(q) c(q[[2]], q[[3]], q[[4]], q[[1]]))
  # Compute the lambdai and Ahat in ascending order.
  qqTs <- lapply(qs, function(q) outer(q, q))
  eig <- eigen(arithmeticMean(qqTs), symmetric=TRUE)
  vals <- rev(eig$values)
  vecs <- cbind(eig$vectors[,4], eig$vectors[,3], eig$vectors[,2], eig$vectors[,1])
  if (vals[[1]] + vals[[4]] <= 0.5) {
    # The data are not bipolar, so we return this garbage p-value function.
    f <- function(r) {NULL}
  }
  else {
    # The data are bipolar. Ensure that vecs has determinant 1.
    if (det(vecs) < 0)
      vecs[,4] <- -vecs[,4]
    # Estimate the center M of the distribution.
    #mHat <- vecs[,4]
    #mHat <- rotationFromQuaternion(c(mHat[[4]], mHat[[1]], mHat[[2]], mHat[[3]]))
    # Approximate the confidence region metric.
    x2s <- lapply(qqTs, function(qqT) diag(t(vecs) %*% qqT %*% vecs))
    cHat <- arithmeticMean(lapply(x2s, function(x2) outer(x2, x2)))
    aHat4By3 <- vecs[c(1, 2, 3, 4), c(1, 2, 3)]
    fHat <- diag(c(
      (vals[[4]] - vals[[1]])^2 / cHat[4, 1],
      (vals[[4]] - vals[[2]])^2 / cHat[4, 2],
      (vals[[4]] - vals[[3]])^2 / cHat[4, 3]))
    fHat <- length(rs) * aHat4By3 %*% fHat %*% t(aHat4By3)
    # Put the metric back into my preferred quaternion order.
    fHat <- rbind(fHat[4,], fHat[1,], fHat[2,], fHat[3,])
    fHat <- cbind(fHat[,4], fHat[,1], fHat[,2], fHat[,3])
    f <- function(r) {
      q <- rotQuaternionFromMatrix(r)
      1 - pchisq(as.numeric(q %*% fHat %*% q), 3)
    }
  }
  f
}



### TANGENT SPACE METHODS ###

#' Sample covariance matrix, approximated in the tangent space at a given rotation.
#'
#' Appropriate only if the sample is tightly concentrated near the center.
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix. Typically the Frechet mean of the rs.
#' @return A 3x3 real matrix (symmetric, non-negative eigenvalues).
rotLeftCovariance <- function(rs, center) {
	vs <- lapply(rs, rotLeftTangentFromMatrix, center)
	ms <- lapply(vs, function(v) {outer(v, v)})
	arithmeticMean(ms)
}

#' Principal geodesic analysis in the tangent space at a given rotation.
#'
#' Appropriate only if the sample is tightly concentrated near the center.
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix. Typically the Frechet mean of the rs.
#' @param numPoints A real number (integer, 0 or >= 3). The number of points to return on each of the six geodesics through the center.
#' @return A list consisting of $magnitudes (3D real vector, nonnegative) and $directions (3x3 real matrix, whose columns are unit-length vectors). The $magnitudes are in decreasing order. The $directions are the corresponding directions, suitable for use in rotMatrixFromLeftTangent. If numPoints >= 1, then there is also a $curves field (list of three lists of (2 numPoints + 1) rotation matrices).
rotLeftPrincipalComponentAnalysis <- function(rs, center, numPoints=0) {
	eig <- eigen(rotLeftCovariance(rs, center), symmetric=TRUE)
  mags <- sqrt(eig$values)
  dirs <- eig$vectors
	result <- list(magnitudes=mags, directions=dirs)
  if (numPoints >= 1) {
    f <- function(s, magDir) {
      rotMatrixFromLeftTangent(s / numPoints * magDir, center)
    }
    curve1 <- lapply(-numPoints:numPoints, f, mags[[1]] * dirs[,1])
    curve2 <- lapply(-numPoints:numPoints, f, mags[[2]] * dirs[,2])
    curve3 <- lapply(-numPoints:numPoints, f, mags[[3]] * dirs[,3])
    result$curves <- list(curve1, curve2, curve3)
  }
  result
}

#' Sampling from the wrapped trivariate normal distribution.
#'
#' @param s A rotation matrix. The center of the distribution.
#' @param kappa A real number (positive). The concentration of the distribution.
#' @param n A real number (positive integer) or NULL.
#' @return Either a rotation matrix (if n is NULL) or a list of n rotation matrices (if n is a number).
rotWrappedTrivariateNormal <- function(s, kappa, n=NULL) {
  if (is.null(n)) {
    v <- rnorm(3, 0, 1 / kappa)
    s %*% rotExp(rotAntisymmetricFromVector(v))
  }
  else
    replicate(n, rotWrappedTrivariateNormal(s, kappa), simplify=FALSE)
}

#' Mahalanobis distance of a rotation relative to a sample.
#'
#' Appropriate only if the sample is tightly clustered and the given rotation is close to its center.
#' @param r A rotation matrix.
#' @param center A rotation matrix. Typically the Frechet mean of a sample.
#' @param covarInv A 3x3 real matrix (symmetric, positive-definite). Typically the inverse of the covariance obtained from rotLeftCovariance.
#' @return A real number.
rotMahalanobisNorm <- function(r, center, leftCovarInv) {
	v <- rotLeftTangentFromMatrix(r, center)
	sqrt(v %*% leftCovarInv %*% v)
}

#' P-value function for any hypothesized mean, using tangent space approximation and Mahalanobis distances.
#' 
#' Appropriate only if the rotations are tightly clustered. May fail if the rotations live on a geodesic curve or surface.
#' @param rs A list of rotation matrices. Typically the result of bootstrapping, MCMC, etc.
#' @param center A rotation matrix. Typically the Frechet mean of the rs.
#' @return A list containing elements $pvalue (an R function from {rotation matrices} to {real numbers}; for any given hypothesized mean, this function produces the p-value for that hypothesis), $center (as passed to this function), $leftCovarInv (real symmetric 3x3 matrix), and $q000, $q025, $q050, $q075, $q095, $q099, $q100 (real numbers; quantiles of Mahalanobis distance). For example, a rotation r with tangent vector v at center is in the 95% confidence region if sqrt(v %*% covarInv %*% v) < q095.
rotMahalanobisPercentiles <- function(rs, center) {
  vs <- lapply(rs, rotLeftTangentFromMatrix, center)
  covar <- arithmeticMean(lapply(vs, function(v) {outer(v, v)}))
  covarInv <- solve(covar)
  norms <- sapply(vs, function(v) {sqrt(v %*% covarInv %*% v)})
  empiricalCDF <- ecdf(norms)
  # Build the p-value function.
  f <- function(r) {
    v <- rotLeftTangentFromMatrix(r, center)
    1 - empiricalCDF(sqrt(v %*% covarInv %*% v))
  }
  # Compute a few popular percentiles.
  qs <- quantile(norms, probs=c(0.00, 0.25, 0.50, 0.75, 0.95, 0.99, 1.00), names=FALSE)
  list(pvalue=f, center=center, leftCovarInv=covarInv,
       q000=qs[[1]], q025=qs[[2]], q050=qs[[3]], q075=qs[[4]], q095=qs[[5]], q099=qs[[6]], q100=qs[[7]])
}

#' Ellipsoidal surface from Mahalanobis inference.
#' 
#' @param center A rotation matrix. Typically the $center from the inference.
#' @param leftCovarInv A 3x3 real matrix (symmetric, positive-definite). Typically the $leftCovarInv from the inference.
#' @param level A real number. Typically $q095^2 from the inference.
#' @param numNonAdapt A real number (non-negative integer). The number of refinements to the sphere that is deformed into the ellipsoid. Incrementing numNonAdapt improves visual quality but increases time and memory requirements by a factor of four.
#' @return A list of triangles, where each triangle is a list of three rotation matrices.
rotEllipsoidTriangles <- function(center, leftCovarInv, level, numNonAdapt=3) {
  # Diagonalize the inverse covariance.
  eig <- eigen(leftCovarInv, symmetric=TRUE)
  q <- eig$vectors
  a <- sqrt(level) * eig$values^(-0.5)
  # Make the ellipsoid in the left-invariant tangent space and transfer it to SO(3).
  sphere <- rayTetrahedralSphere(numNonAdapt)
  ellipsoid <- lapply(sphere, function(tri) lapply(tri, function(v) rotMatrixFromLeftTangent(q %*% (a * v), center)))
  ellipsoid
}

#' Inference about the population mean, based on non-parametric bootstrapping.
#' 
#' This function bootstraps the rotation mean, returning two pieces of information. The first is a list of the bootstrapped means. The user should rotEqualAnglePlot or rotEqualVolumePlot this list, to make sure that the means form a fairly tight ellipsoidal cluster. If so, then the second piece of information may be used: An R function that takes a rotation R0 as input, as produces as output a p-value for the hypothesis that the mean of the population is R0. This p-value is the fraction of means that are farther from the mean of means than R0 is, based on Mahalanobis distance in the set of means.
#' @param rs A list of rotation matrices.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @param func An R function. Either rotFrechetMean or rotProjectedMean (or something equivalent).
#' @return A list consisting of $bootstraps (a list of rotation matrices) and everything returned by rotMahalanobisPercentiles.
rotBootstrapInference <- function(rs, numBoots, func=rotFrechetMean) {
  boots <- replicate(numBoots, func(sample(rs, length(rs), replace=TRUE)), simplify=FALSE)
  bootMean <- func(boots)
  infer <- rotMahalanobisPercentiles(boots, bootMean)
  infer$bootstraps <- boots
  infer
}

#' Two-sample inference about the difference in population means, based on non-parametric bootstrapping.
#' 
#' This function bootstraps the difference of the means of the two samples, returning two pieces of information. The first is a list of the bootstrapped differences. The user should rotEqualAnglePlot or rotEqualVolumePlot this list, to make sure that the differences form a fairly tight ellipsoidal cluster. If so, then the second piece of information may be used: An R function that takes a rotation R0 as input, as produces as output a p-value for the hypothesis that the difference of means is R0. This p-value is the fraction of differences that are farther from the mean of differences than R0 is, based on the Mahalanobis distance of the differences.
#' @param firsts A list of rotation matrices.
#' @param seconds A list of rotation matrices.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @param func An R function. Either rotFrechetMean or rotProjectedMean (or something equivalent).
#' @return A list consisting of $bootstraps (a list of rotation matrices) and everything returned by rotMahalanobisPercentiles.
rotTwoSampleBootstrapInference <- function(firsts, seconds, numBoots, func=rotFrechetMean) {
  f <- function() {
    firstMean <- func(sample(firsts, length(firsts), replace=TRUE))
    secondMean <- func(sample(seconds, length(seconds), replace=TRUE))
    secondMean %*% t(firstMean)
  }
  boots <- replicate(numBoots, f(), simplify=FALSE)
  bootMean <- func(boots, ...)
  infer <- rotMahalanobisPercentiles(boots, bootMean)
  infer$bootstraps <- boots
  infer
}



# Simulation arms for the accompanying simulation studies, and the stage-two
# diagnostic. Not part of the MRKCFC method; kept in the package so the
# version that produces the reported results is the version that ships.
#
# Six arms ordered by how many moments exist, which is the axis Proposition 1
# of the manuscript turns on. The Fleishman arms are read as (skewness, excess
# kurtosis): skewness 4 would force an excess kurtosis of at least 14.
#
#   1  Gaussian                        skew 0    ex.kurt 0       all moments
#   2  Fleishman(0, 4)                 skew 0    ex.kurt 4       all moments
#   3  Fleishman(1.5, 4)               skew 1.5  ex.kurt 4       all moments
#   4  Contaminated N, eps=.05, k=5    skew 0    ex.kurt 16.96   all moments
#   5  t, nu = 3                       skew 0    ex.kurt inf     variance finite
#   6  t, nu = 1.5                     skew 0    ex.kurt inf     variance infinite
#
# Every arm is calibrated so each score component has the same IQR as the
# Gaussian arm. Arms are matched on a robust dispersion measure, not on
# variance, because no common variance exists across the design. Without this
# a reader can attribute any difference to signal-to-noise rather than tails.
#
# The constants below came from numerical solves checked against simulation
# and re-verified analytically: the Fleishman coefficients solve the moment
# equations to 3e-5, and IQR_SCALE matches exact quantile-based values. Do not
# recompute or round them.

FLEISHMAN <- list(
  `2` = c(C0 = 0.0,       C1 = 0.737381, C2 = 0.0,       C3 = 0.080925),
  `3` = c(C0 = -0.211594, C1 = 0.846763, C2 = 0.211594,  C3 = 0.034551)
)

# IQR of each standardized arm relative to a unit-variance normal.
IQR_SCALE <- c(`1` = 1.0000, `2` = 1.2909, `3` = 1.1584,
               `4` = 0.9550, `5` = 0.8818, `6` = 0.7730)

#' @title Vale and Maurelli intermediate correlation matrix
#'
#' @description
#' The Fleishman polynomial applied marginally to correlated normal variates
#' attenuates the correlations. This solves for the intermediate normal
#' correlation that yields the target correlation after transformation.
#'
#' @details
#' The correction is small when the target correlations are weak, but it grows
#' with both the correlation and the number of components, so the equation is
#' solved rather than tabulated.
#'
#' @param R The target correlation matrix.
#' @param arm Either \code{2} or \code{3}, naming the Fleishman arm.
#'
#' @return A correlation matrix of the same dimension as \code{R}, to be passed
#'   to the normal generator.
#'
#' @references
#' Vale, C. D. and Maurelli, V. A. (1983). Simulating multivariate nonnormal
#' distributions. \emph{Psychometrika}, \bold{48}(3), 465-471.
#'
#' @seealso \code{\link{gen_scores}}
#'
#' @export
vm_intermediate <- function(R, arm) {
  k <- FLEISHMAN[[as.character(arm)]]
  b <- k[["C1"]]; cc <- k[["C2"]]; d <- k[["C3"]]
  f <- function(rs, rt) rs * (b^2 + 6*b*d + 9*d^2) + rs^2 * (2*cc^2) + rs^3 * (6*d^2) - rt
  out <- diag(nrow(R))
  for (i in seq_len(nrow(R) - 1L)) for (j in (i + 1L):nrow(R)) {
    out[i, j] <- out[j, i] <-
      if (R[i, j] == 0) 0 else uniroot(f, c(-0.999, 0.999), rt = R[i, j])$root
  }
  out
}

rmvn <- function(n, S) matrix(rnorm(n * nrow(S)), n) %*% chol(S)

#' @title Generate functional principal component scores for one simulation arm
#'
#' @description
#' Draws \emph{n} score vectors from one of six distributions ordered by how many
#' moments exist, each calibrated to the interquartile range of the Gaussian arm
#' so that tail behaviour is not confounded with signal-to-noise.
#'
#' @details
#' Arm 1 is Gaussian. Arms 2 and 3 apply a Fleishman power-method polynomial to
#' correlated normal variates, with the intermediate correlation from
#' \code{\link{vm_intermediate}}, giving (skewness, excess kurtosis) of (0, 4)
#' and (1.5, 4). Neither is a skew-normal variate: the skew-normal family attains
#' at most an excess kurtosis of about 0.869.
#'
#' Arm 4 is a scale-contaminated normal: a proportion \code{eps} of subjects have
#' their scores inflated by a factor \code{k}. Contamination is by scale rather
#' than by a mean shift, because a shifted contaminant has no well-defined true
#' class label and would confound the clustering metrics with the definition of
#' an outlier.
#'
#' Arms 5 and 6 are multivariate \emph{t} with 3 and 1.5 degrees of freedom. In
#' arm 6 no covariance exists, so \code{Z} is the \emph{scatter} matrix there
#' rather than the covariance.
#'
#' @param n Number of subjects.
#' @param arm An integer from 1 to 6; see Details.
#' @param Z The score dispersion matrix: the covariance in arms 1 to 5, the
#'   scatter matrix in arm 6.
#' @param eps Contamination proportion, arm 4 only.
#' @param k Contamination scale inflation, arm 4 only.
#'
#' @return An \emph{n} by \code{ncol(Z)} matrix of scores.
#'
#' @references
#' Fleishman, A. I. (1978). A method for simulating non-normal distributions.
#' \emph{Psychometrika}, \bold{43}(4), 521-532.
#'
#' @examples
#' set.seed(1)
#' Z <- diag(c(9, 4))
#' s <- gen_scores(2000, arm = 5, Z = Z)
#' # IQR matching: comparable spread to the Gaussian arm despite infinite kurtosis
#' rbind(t3 = apply(s, 2, stats::IQR),
#'       gaussian = apply(gen_scores(2000, arm = 1, Z = Z), 2, stats::IQR))
#'
#' @export
gen_scores <- function(n, arm, Z, eps = 0.05, k = 5) {
  sd <- sqrt(diag(Z))
  R  <- Z / outer(sd, sd)

  X <- switch(as.character(arm),
    "1" = rmvn(n, R),
    "2" = ,
    "3" = {
      z <- rmvn(n, vm_intermediate(R, arm))
      cf <- FLEISHMAN[[as.character(arm)]]
      cf[["C0"]] + cf[["C1"]] * z + cf[["C2"]] * z^2 + cf[["C3"]] * z^3
    },
    "4" = {
      x <- rmvn(n, R)
      sel <- runif(n) < eps          # one mask, drawn once
      x[sel, ] <- x[sel, ] * k
      x
    },
    "5" = ,
    "6" = {
      nu <- if (arm == 5) 3 else 1.5
      rmvn(n, R) * sqrt(nu / rchisq(n, nu))    # unit scatter; IQR_SCALE calibrates
    },
    stop("arm must be 1..6")
  )
  X * IQR_SCALE[[as.character(arm)]] * rep(sd, each = n)
}


#' @title Principal angles between two subspaces
#'
#' @description
#' The angles, in degrees, between the column spans of two matrices. Used to
#' compare an estimated eigenfunction subspace against the truth, or two
#' estimators against each other, without depending on the sign or the ordering
#' of individual eigenvectors.
#'
#' @param A A numeric matrix whose columns span the first subspace.
#' @param B A numeric matrix whose columns span the second subspace.
#'
#' @return A numeric vector of angles in degrees, in increasing order.
#'
#' @examples
#' A <- diag(3)[, 1:2]
#' B <- cbind(c(1, 0, 0), c(0, cos(pi / 6), sin(pi / 6)))
#' principal_angles(A, B)
#'
#' @export
principal_angles <- function(A, B) {
  Qa <- qr.Q(qr(A)); Qb <- qr.Q(qr(B))
  acos(pmin(pmax(svd(crossprod(Qa, Qb))$d, -1), 1)) * 180 / pi
}

#' @title Angle between the moment and sign second-stage subspaces
#'
#' @description
#' The largest principal angle between the leading \code{M_keep} eigenvectors
#' of the second-moment matrix of a score matrix and those of its spatial-sign
#' covariance matrix. Near zero means the two second stages agree on the
#' subspace; a large value means the sign scatter sees a different structure
#' from the second moments, as happens under contamination or heavy tails.
#'
#' @param Xi A numeric \emph{n} by \emph{M} matrix of stacked principal component
#'   scores.
#' @param M_keep Number of leading components to compare.
#'
#' @return A single angle in degrees.
#'
#' @seealso \code{\link{robust_stage2}}, \code{\link{principal_angles}}
#'
#' @examples
#' set.seed(1)
#' Z <- diag(c(9, 4, 1))
#' stage2_diagnostic(gen_scores(500, arm = 1, Z = Z), M_keep = 2)
#' stage2_diagnostic(gen_scores(500, arm = 6, Z = Z), M_keep = 2)
#'
#' @export
stage2_diagnostic <- function(Xi, M_keep = 3L) {
  max(principal_angles(robust_stage2(Xi, M_keep, robust = FALSE)$vectors,
                       robust_stage2(Xi, M_keep, robust = TRUE)$vectors))
}

# Robust second stage for MFPCA, replacing the sample second-moment
# eigendecomposition of the stacked score matrix.
#
# The original algorithm forms Z-hat = n^-1 t(Xi) %*% Xi and takes its
# eigenvectors. That is a sample second-moment matrix: unbounded influence
# function, so robustness attained in the marginal Kendall FPCA step is lost
# here: the method was robust in its marginal step only.
#
# Replacement: the spatial sign covariance matrix of the centred score
# vectors. Under an elliptical law its eigenvectors are consistent for those
# of the scatter matrix (Visuri, Koivunen & Oja 2000; Magyar & Tyler 2014),
# which is the same argument Zhong et al. use marginally. Eigenvalues are NOT
# preserved, so scales are recovered separately by a robust scale of the
# scores projected onto each estimated direction.
#
# The pairwise multivariate Kendall's tau of Han & Liu (2018) shares the
# eigenvector property but costs O(n^2 p^2) per call. Inside the K-centers
# loop that is prohibitive; the spatial-sign form is O(n p^2).

#' @title Spatial median of a score matrix
#'
#' @description
#' The multivariate \eqn{L_1} median, computed by Weiszfeld's iteration from a
#' coordinatewise-median start. Unlike the mean it has a population target
#' without any moment assumption, which is what makes it usable on heavy-tailed
#' principal component scores.
#'
#' @param X A numeric matrix, one row per subject.
#' @param tol Convergence tolerance on the change in the estimate.
#' @param max_it Maximum number of iterations.
#'
#' @return A numeric vector of length \code{ncol(X)}.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(stats::rt(200, df = 2), ncol = 2)
#' spatial_median(X)
#'
#' @export
spatial_median <- function(X, tol = 1e-8, max_it = 200L) {
  m <- apply(X, 2, median)
  for (i in seq_len(max_it)) {
    d <- sqrt(rowSums((X - rep(m, each = nrow(X)))^2))
    d[d < 1e-12] <- 1e-12
    m_new <- colSums(X / d) / sum(1 / d)
    if (sqrt(sum((m_new - m)^2)) < tol) return(m_new)
    m <- m_new
  }
  m
}

#' @title Spatial sign covariance matrix
#'
#' @description
#' The covariance of the unit vectors obtained by projecting deviations from the
#' spatial median onto the unit sphere. Under an elliptical law its eigenvectors
#' are consistent for those of the scatter matrix, which is the same argument
#' Kendall's tau uses at the operator level, so no moment assumption is required.
#'
#' @details
#' Eigenvectors are recovered; eigenvalues are not. The matrix has trace 1 by
#' construction, so its eigenvalues carry no information about the spread along
#' each direction. \code{\link{robust_stage2}} recovers that separately.
#'
#' @param X A numeric matrix, one row per subject.
#'
#' @return A \code{ncol(X)} by \code{ncol(X)} symmetric matrix.
#'
#' @references
#' Visuri, S., Koivunen, V. and Oja, H. (2000). Sign and rank covariance
#' matrices. \emph{Journal of Statistical Planning and Inference}, \bold{91}(2),
#' 557-575.
#'
#' @seealso \code{\link{robust_stage2}}
#'
#' @export
sscm <- function(X) {
  U <- X - rep(spatial_median(X), each = nrow(X))
  nrm <- sqrt(rowSums(U^2))
  keep <- nrm > 1e-12
  U <- U[keep, , drop = FALSE] / nrm[keep]
  crossprod(U) / nrow(U)
}

#' @title Second-stage eigenstructure of the stacked scores
#'
#' @description
#' The rotation of the second stage of the multivariate FPCA: eigenvectors of
#' the spatial sign covariance matrix of the stacked univariate scores, with the
#' spread along each estimated direction recovered by a robust scale rather than
#' by a second moment. With \code{robust = FALSE} it reduces to the ordinary
#' sample second-moment eigendecomposition, which is the comparator.
#'
#' @details
#' The spatial sign covariance matrix has trace 1, so inverting or ranking by its
#' eigenvalues would put the result on an arbitrary scale. The spread along each
#' estimated direction is therefore recovered as \eqn{(\mathrm{IQR}/1.349)^2} of
#' the projected scores, which is Fisher-consistent for the variance at the
#' normal and exists without a moment assumption.
#'
#' The components are ordered by that recovered scale, not by the sign
#' covariance eigenvalue. In the population the two orderings agree, since the
#' sign eigenvalues are monotone in the scatter eigenvalues, but in a finite
#' sample under heavy tails they can differ, and the reconstruction truncates to
#' the leading components, so the order decides which subspace is used.
#'
#' @param Xi A numeric \emph{n} by \emph{M} matrix of stacked principal
#'   component scores.
#' @param M_keep Number of components to retain; defaults to all.
#' @param robust Logical; \code{TRUE} uses the spatial sign covariance matrix,
#'   \code{FALSE} the sample second-moment matrix.
#'
#' @return A list with components
#'   \item{vectors}{an \emph{M} by \code{M_keep} matrix of eigenvectors in columns}
#'   \item{values}{the spread along each retained direction, on the score scale}
#'
#' @seealso \code{\link{sscm}}, \code{\link{whiten_scores}},
#'   \code{\link{MFPCA_k}}
#'
#' @examples
#' set.seed(1)
#' Xi <- matrix(stats::rt(300, df = 2), ncol = 3) %*% diag(c(3, 2, 1))
#' robust_stage2(Xi, M_keep = 2)$values
#'
#' @export
robust_stage2 <- function(Xi, M_keep = ncol(Xi), robust = TRUE) {
  stopifnot(is.matrix(Xi), M_keep >= 1L, M_keep <= ncol(Xi))
  S <- if (robust) sscm(Xi) else crossprod(Xi) / nrow(Xi)
  e <- eigen(S, symmetric = TRUE)

  if (robust) {
    # SSCM eigenvalues are not the score variances. Recover scale along each
    # estimated direction with a robust scale. IQR/1.349 is Fisher-consistent
    # for the standard deviation at the normal, and matches the estimator now
    # used in robust_normalize().
    #
    # ORDER BY THAT RECOVERED SCALE, not by the SSCM eigenvalue. In the
    # population the two orderings agree, since the SSCM eigenvalues are
    # monotone in the scatter eigenvalues, but in a finite sample under heavy
    # tails they do not. On simulation arm 6, which has infinite variance, the
    # SSCM order returned spreads of 8.54, 8.78, 7.87, ... : component 2 was
    # more important than component 1.
    #
    # This matters because Step 3(d) truncates to the LEADING components, so a
    # wrong order silently reconstructs from the wrong subspace. A comparison
    # of principal angles between spans cannot see it, being order-invariant
    # by construction; it shows only once the reconstruction truncates.
    nu_all <- apply(Xi %*% e$vectors, 2, function(z) (IQR(z) / 1.349)^2)
    keep   <- order(nu_all, decreasing = TRUE)[seq_len(M_keep)]
    C  <- e$vectors[, keep, drop = FALSE]
    nu <- nu_all[keep]
  } else {
    # eigen() already returns decreasing eigenvalues, so the moment branch is
    # ordered by construction.
    C  <- e$vectors[, seq_len(M_keep), drop = FALSE]
    nu <- e$values[seq_len(M_keep)]
  }
  list(vectors = C, values = nu)
}


# Single entry point for the stage-two rotation.
#
# Every site that needs the rotation calls this: MFPCA_k(), and both branches
# of predicted_yhat(). Computing it in more than one place invites the two
# sites to disagree on the sign of a component, since eigenvectors are defined
# only up to sign and different solvers choose differently; one function
# removes that possibility, and .mrkcfc_mfpca_reconstruct() below is in any
# case invariant to the sign.
#
# mfpca_scatter names the scatter estimator, not a vague robustness property.
# "scatter" rather than "covariance" is deliberate: under the heaviest-tailed
# simulation arm no covariance exists, and that distinction is what
# Proposition 1 turns on.
#
# Both branches use eigen(S, symmetric = TRUE). Z-hat is symmetric by
# construction, Step 3(c) specifies an eigendecomposition and needs the
# eigenvalues nu-hat, and the symmetric solver guarantees real eigenvalues.
# The moment branch divides by n, as the manuscript defines Z-hat.

.mrkcfc_mfpca_rotation <- function(full_xiEst, mfpca_scatter = c("sign", "moment")) {
  mfpca_scatter <- match.arg(mfpca_scatter)
  robust_stage2(full_xiEst, robust = identical(mfpca_scatter, "sign"))$vectors
}


# Step 3(d) and Equation (14): multivariate eigenfunctions, multivariate
# scores, and the class-specific prediction for one subject.
#
# Follows Happ and Greven (2018), JASA 113(522), Section 3.2 step 4:
#
#   psi_m^(j)(t) = sum_{n=1}^{M_j} [c_m]_n^(j) phi_n^(j)(t),  m = 1,...,M+   (9)
#   rho_{i,m}    = sum_{j=1}^{p} sum_{n=1}^{M_j} [c_m]_n^(j) xi_{i,n}^(j)
#                = Xi_{i,.} c_m                                             (10)
#
# The score is ONE scalar per component, formed from the subject's FULL
# stacked score vector across every element. That cross-element coupling is
# the entire point of the second stage; a per-element score is rotated
# univariate FPCA.
#
# Truncation follows Equation (14): element j is reconstructed from the
# LEADING M_j multivariate components, not from components index_j. Retaining
# all M+ would be lossless and pointless, since C_j %*% t(C) is then just the
# index_j rows of the identity and the reconstruction collapses back to the
# univariate one.
#
# This form is invariant to eigenvector sign by construction: under C -> C D
# both rotation[index, keep] and rho[keep] pick up D, and D %*% D = I.

.mrkcfc_mfpca_reconstruct <- function(rotation, pi_obs, xi_i, n_M, mu_obs) {
  J <- length(pi_obs)
  full_xi_i <- unlist(lapply(xi_i, as.numeric), use.names = FALSE)
  stopifnot(length(full_xi_i) == nrow(rotation))

  rho_i <- as.vector(crossprod(rotation, full_xi_i))          # Equation (10)

  fittedX <- multi_eigen <- vector("list", J)
  index <- 0L
  for (j in seq_len(J)) {
    index <- seq.int(from = max(index) + 1L, length.out = n_M[[j]])
    keep  <- seq_len(n_M[[j]])
    Psi_j <- pi_obs[[j]] %*% rotation[index, keep, drop = FALSE]  # Equation (9)
    multi_eigen[[j]] <- t(Psi_j)
    fittedX[[j]] <- as.numeric(mu_obs[[j]] + Psi_j %*% rho_i[keep])
  }
  # multi_score is a single vector of length M+, not one per element.
  list(fittedX = fittedX, multi_score = rho_i, multi_eigen = multi_eigen)
}

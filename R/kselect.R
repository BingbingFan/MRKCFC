# Cluster-validity criteria for selecting K.
#
# WHY THESE WORK ON SCORES AND NOT ON CURVES
#
# A Mahalanobis distance needs a common coordinate system in which to form a
# scatter matrix, and sparse irregularly observed curves do not have one:
# subject i is seen at its own times and subject j at different ones, so the
# raw observations are not comparable coordinates. The multivariate FPC scores
# are exactly the missing coordinate system -- M+ numbers per subject, the same
# M+ for everyone, and the space the K-centers loop actually clusters in. This
# is a constraint, not a preference.
#
# WHY THE WHITENING IS ESTIMATED ONCE, GLOBALLY
#
# The obvious construction -- whiten by the pooled WITHIN-cluster scatter of
# each candidate partition -- is degenerate, and silently so. If S is the
# pooled within-cluster scatter of the same partition, then the within-cluster
# sum of squares in the whitened metric is
#
#     sum_k sum_{i in k} d2_S(x_i, m_k) = trace(S^-1 W) = trace(n I) = n p
#
# identically, for every K. W would be constant, CH would reduce to B alone,
# and B alone is monotone in K, so the criterion would always choose the
# largest K offered. The metric must therefore be fixed across the K it is
# comparing. Here it is estimated once from the unpartitioned scores, so
# "mahalanobis" is a whitening applied before the criteria rather than a term
# inside them.
#
# WHY IT IS ROBUST
#
# With robust = TRUE both the location and the scatter are computed without
# moment assumptions -- spatial median for location, C diag(nu) t(C) from
# robust_stage2() for scatter -- so the criteria stay defined on simulation
# arms 5 and 6, whose scores have infinite fourth moment and infinite variance
# respectively.

#' @title Whiten scores to a Mahalanobis metric
#'
#' @description
#' Returns coordinates in which Euclidean distance equals Mahalanobis distance in
#' the original score space, using a robust location and scatter by default.
#'
#' @details
#' \eqn{S^{-1/2} = C \,\mathrm{diag}(\nu^{-1/2})\, C^{T}}, but the trailing
#' \eqn{C^{T}} is an orthogonal rotation and orthogonal maps preserve Euclidean
#' distance, so it is dropped: the returned coordinates differ from the textbook
#' ones by a rotation and give identical distances.
#'
#' A direction with near-zero spread carries no information and would be
#' amplified without bound by \eqn{\nu^{-1/2}}, so it is floored against the
#' largest and thereby neutralised rather than allowed to dominate.
#'
#' @param scores A numeric \emph{n} by \emph{M} matrix of principal component
#'   scores.
#' @param robust Logical; robust location and scatter, or their moment analogues.
#' @param tol Relative floor below which a direction's spread is treated as
#'   degenerate.
#'
#' @return A numeric matrix with the same dimensions as \code{scores}.
#'
#' @seealso \code{\link{robust_stage2}}, \code{\link{select_K}}
#'
#' @export
whiten_scores <- function(scores, robust = TRUE, tol = 1e-8) {
  stopifnot(is.matrix(scores), nrow(scores) >= 2L)
  centre <- if (robust) spatial_median(scores) else colMeans(scores)
  Xc <- scores - rep(centre, each = nrow(scores))
  e <- robust_stage2(Xc, robust = robust)

  # A direction with (near) zero spread carries no information and would be
  # amplified without bound by nu^-1/2. Floor it against the largest, so a
  # degenerate component is neutralised rather than allowed to dominate.
  nu <- e$values
  nu[!is.finite(nu) | nu < tol * max(nu, na.rm = TRUE)] <- tol * max(nu, na.rm = TRUE)
  Xc %*% e$vectors %*% diag(1 / sqrt(nu), nrow = length(nu))
}

# Dispatch on the distance argument once, so both criteria below see a plain
# matrix in which Euclidean distance is the intended distance.
.mrkcfc_kselect_space <- function(scores, distance = c("euclidean", "mahalanobis"),
                                 robust = TRUE) {
  distance <- match.arg(distance)
  if (distance == "euclidean") as.matrix(scores)
  else whiten_scores(as.matrix(scores), robust = robust)
}

# Calinski-Harabasz index.
#
#     CH(K) = [B(K)/(K-1)] / [W(K)/(n-K)]
#
# with B(K) = sum_k n_k ||m_k - m||^2 and W(K) = sum_k sum_{i in k} ||x_i - m_k||^2,
# exactly as Calinski and Harabasz (1974) define it.
#
# THE DEGREES-OF-FREEDOM DIVISORS ARE NOT COSMETIC. W falls and B rises with
# every extra cluster, so the bare ratio B/W increases monotonically in K and
# its argmax is always the largest K offered. (K-1) and (n-K) are what turn a
# monotone quantity into one with an interior maximum, which is the entire
# basis for using it to select K; tests/testthat/test-kselect.R demonstrates
# the difference.
#
# `robust` chooses the cluster locations: spatial medians rather than means,
# which matters on the heavy-tailed arms where a mean has no population target.

#' @title Calinski-Harabasz index on principal component scores
#'
#' @description
#' \deqn{CH(K) = \frac{B(K)/(K-1)}{W(K)/(n-K)}}
#' with \eqn{B(K) = \sum_k n_k \|m_k - m\|^2} and
#' \eqn{W(K) = \sum_k \sum_{i \in k} \|x_i - m_k\|^2}, as Calinski and Harabasz
#' (1974) define it.
#'
#' @details
#' The degrees-of-freedom divisors are not cosmetic. \eqn{W} falls and \eqn{B}
#' rises with every extra cluster, so the bare ratio \eqn{B/W} increases
#' monotonically in \eqn{K} and its argmax is always the largest \eqn{K} offered.
#' \eqn{(K-1)} and \eqn{(n-K)} are what turn a monotone quantity into one with an
#' interior maximum, which is the entire basis for using it to select \eqn{K}.
#'
#' \code{robust} chooses the cluster locations: spatial medians rather than
#' means, which matters where a mean has no population target.
#'
#' @param scores A numeric \emph{n} by \emph{M} matrix of principal component
#'   scores from one global fit.
#' @param labels A partition of the \emph{n} subjects.
#' @param distance \code{"euclidean"} on the raw scores, or \code{"mahalanobis"}
#'   after \code{\link{whiten_scores}}.
#' @param robust Logical; robust locations and scatter.
#'
#' @return The index, or \code{NA} when it is undefined (fewer than two clusters,
#'   as many clusters as subjects, or zero within-cluster spread).
#'
#' @references
#' Calinski, T. and Harabasz, J. (1974). A dendrite method for cluster analysis.
#' \emph{Communications in Statistics}, \bold{3}(1), 1-27.
#'
#' @seealso \code{\link{select_K}}, \code{\link{silhouette_index}}
#'
#' @export
ch_index <- function(scores, labels, distance = c("euclidean", "mahalanobis"),
                     robust = TRUE) {
  X <- .mrkcfc_kselect_space(scores, distance, robust)
  labels <- as.integer(as.factor(labels))
  cls <- sort(unique(labels))
  K <- length(cls); n <- nrow(X)
  if (K < 2L || K >= n) return(NA_real_)

  loc <- function(M) if (robust && nrow(M) > 1L) spatial_median(M) else colMeans(M)
  m <- loc(X)
  B <- 0; W <- 0
  for (k in cls) {
    Xk <- X[labels == k, , drop = FALSE]
    mk <- loc(Xk)
    B <- B + nrow(Xk) * sum((mk - m)^2)
    W <- W + sum((Xk - rep(mk, each = nrow(Xk)))^2)
  }
  if (W <= 0) return(NA_real_)
  (B / (K - 1L)) / (W / (n - K))
}

# Silhouette index.
#
# Rousseeuw (1987):
#
#     a(i) = mean distance from i to the OTHER members of its own cluster
#     b(i) = min over other clusters of the mean distance from i to that cluster
#     s(i) = (b(i) - a(i)) / max(a(i), b(i)),  s(i) = 0 for a singleton cluster
#
# method = "pairwise" is that definition. method = "centroid" is the simplified
# silhouette, which replaces both means by the distance to the cluster location;
# it is far cheaper and is what silhouette_score() computes on curves, where
# pairwise distances between irregularly observed subjects are not directly
# available. The two can disagree in sign for boundary subjects, so the choice
# is recorded rather than assumed.
#
# Note that this is a DISTANCE, not a squared distance. silhouette_score()
# scores on squared L2 deviations, and s() is not invariant to that: squaring
# exaggerates large gaps, so (b^2-a^2)/max(a^2,b^2) reads as more separated
# than (b-a)/max(a,b) on the same data.

#' @title Silhouette index on principal component scores
#'
#' @description
#' Rousseeuw's silhouette computed in score space, as an alternative to the
#' curve-space \code{\link{silhouette_score}}.
#'
#' @details
#' With \code{method = "pairwise"} this is Rousseeuw's (1987) definition:
#' \eqn{a(i)} is the mean distance from \emph{i} to the other members of its own
#' cluster, \eqn{b(i)} the minimum over other clusters of the mean distance from
#' \emph{i} to that cluster, and \eqn{s(i) = (b_i - a_i)/\max(a_i, b_i)}, with
#' \eqn{s(i) = 0} for a singleton cluster.
#'
#' With \code{method = "centroid"} it is the simplified silhouette, replacing
#' both means by the distance to the cluster location. That is far cheaper and is
#' the form \code{\link{silhouette_score}} computes on curves, where pairwise
#' distances between irregularly observed subjects are not directly available.
#' The two can disagree in sign for boundary subjects, so the choice is recorded
#' rather than assumed.
#'
#' Note that this is a distance, not a squared distance, and the silhouette is
#' not invariant to that. \code{\link{silhouette_score}} scores on squared
#' \eqn{L_2} deviations; squaring exaggerates large gaps, so
#' \eqn{(b^2-a^2)/\max(a^2,b^2)} reads as more separated than
#' \eqn{(b-a)/\max(a,b)} on the same data. The two are not interchangeable.
#'
#' @param scores A numeric \emph{n} by \emph{M} matrix of principal component
#'   scores from one global fit.
#' @param labels A partition of the \emph{n} subjects.
#' @param distance \code{"euclidean"} or \code{"mahalanobis"}.
#' @param robust Logical; robust locations and scatter.
#' @param method \code{"pairwise"} for Rousseeuw's definition, \code{"centroid"}
#'   for the simplified silhouette.
#'
#' @return The mean silhouette over subjects, or \code{NA} with fewer than two
#'   clusters.
#'
#' @references
#' Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the interpretation
#' and validation of cluster analysis. \emph{Journal of Computational and Applied
#' Mathematics}, \bold{20}, 53-65.
#'
#' @seealso \code{\link{select_K}}, \code{\link{silhouette_score}}
#'
#' @export
silhouette_index <- function(scores, labels, distance = c("euclidean", "mahalanobis"),
                             robust = TRUE, method = c("pairwise", "centroid")) {
  method <- match.arg(method)
  X <- .mrkcfc_kselect_space(scores, distance, robust)
  labels <- as.integer(as.factor(labels))
  cls <- sort(unique(labels))
  K <- length(cls); n <- nrow(X)
  if (K < 2L) return(NA_real_)

  if (method == "centroid") {
    loc <- function(M) if (robust && nrow(M) > 1L) spatial_median(M) else colMeans(M)
    cent <- lapply(cls, function(k) loc(X[labels == k, , drop = FALSE]))
    D <- vapply(cent, function(mk) sqrt(rowSums((X - rep(mk, each = n))^2)), numeric(n))
  } else {
    D0 <- as.matrix(stats::dist(X))
    # Column k holds each subject's mean distance to cluster k. For a subject's
    # OWN cluster the diagonal zero must be excluded and the divisor is n_k - 1,
    # which is what separates a(i) from a plain cluster mean.
    D <- vapply(cls, function(k) {
      ink <- labels == k
      s <- rowSums(D0[, ink, drop = FALSE])
      nk <- sum(ink)
      ifelse(ink, if (nk > 1L) s / (nk - 1L) else NA_real_, s / nk)
    }, numeric(n))
  }

  own <- cbind(seq_len(n), match(labels, cls))
  a <- D[own]
  b <- vapply(seq_len(n), function(i) min(D[i, -own[i, 2L]]), numeric(1))
  s <- (b - a) / pmax(a, b)
  s[!is.finite(s)] <- 0          # singleton cluster: s(i) = 0 by definition
  mean(s)
}

# Sweep K and return both criteria.
#
# `scores` must come from ONE global fit shared across the K being compared. A
# criterion compares partitions of a fixed geometry; recomputing the scores per
# K would move the space at the same time as the partition and the comparison
# would no longer be about K.

#' @title Sweep the number of classes and return both criteria
#'
#' @description
#' Evaluates \code{\link{ch_index}} and \code{\link{silhouette_index}} for a list
#' of candidate partitions, so that agreement or disagreement between two
#' criteria is visible rather than a single weak criterion carrying the decision.
#'
#' @details
#' \code{scores} must come from \emph{one} global fit shared across the \code{K}
#' being compared, for instance \code{\link{MFPCA_k}} with a single class over
#' all subjects. A criterion compares partitions of a fixed geometry;
#' recomputing the scores per \code{K} would move the space at the same time as
#' the partition, and the comparison would no longer be about \code{K}.
#'
#' @param scores A numeric \emph{n} by \emph{M} matrix of multivariate principal
#'   component scores from one global fit.
#' @param partitions A list of candidate partitions, typically the
#'   \code{res_clust$g} columns from \code{\link{MKCFC}} fits over a range of
#'   \code{K}.
#' @param distance \code{"euclidean"} or \code{"mahalanobis"}.
#' @param robust Logical; robust locations and scatter.
#' @param silhouette_method Passed to \code{\link{silhouette_index}}.
#'
#' @return A data frame with one row per partition and columns \code{K},
#'   \code{CH}, \code{silhouette}, \code{distance} and \code{robust}.
#'
#' @seealso \code{\link{silhouette_score}} for the curve-space criterion.
#'
#' @examples
#' set.seed(1)
#' scores <- rbind(matrix(stats::rnorm(200, -2), ncol = 2),
#'                 matrix(stats::rnorm(200,  2), ncol = 2))
#' parts <- lapply(2:4, function(k) stats::kmeans(scores, centers = k)$cluster)
#' select_K(scores, parts, robust = FALSE)
#'
#' @export
select_K <- function(scores, partitions, distance = c("euclidean", "mahalanobis"),
                     robust = TRUE, silhouette_method = c("pairwise", "centroid")) {
  distance <- match.arg(distance)
  silhouette_method <- match.arg(silhouette_method)
  Ks <- vapply(partitions, function(g) length(unique(g)), integer(1))
  data.frame(
    K = Ks,
    CH = vapply(partitions, ch_index, numeric(1),
                scores = scores, distance = distance, robust = robust),
    silhouette = vapply(partitions, silhouette_index, numeric(1),
                        scores = scores, distance = distance, robust = robust,
                        method = silhouette_method),
    distance = distance, robust = robust, stringsAsFactors = FALSE
  )
}

#
# Step 2 of the algorithm:
#
#   "Initialization: Perform MRFPCA on all subjects to obtain the marginal mean
#    function mu(t), eigenfunctions psi(t) and multivariate FPC scores xi_i.
#    Cluster these multivariate FPC scores into K groups using K-means or
#    randomly assign the class memberships ... Similarly, perform MFPCA when
#    estimating the MKCFC model."
#
# Why the start follows the same estimator as the loop. The K-centers loop is
# a Lloyd-type alternating scheme that moves to a LOCAL optimum, so the
# starting partition selects the basin and the basin selects the answer. The
# initialiser is therefore a component of the estimator, and the breakdown
# point of a composite estimator is bounded above by the minimum over its
# components: a non-robust start means a single extreme observation can
# determine the partition no matter how robust every later stage is. Under a
# score distribution without finite fourth moments the empirical covariance
# surface of a PACE-type FPCA is not even consistent, and under infinite
# variance it has no population target.
#
# `mfpca_scatter` is threaded through so the initialiser's second stage matches
# the loop's. Otherwise a robust loop would start from a moment-based rotation
# and reintroduce the weakest-link problem in miniature.
#
# Note on k-means itself: nstart is left at the R default of 1, as the
# algorithm specifies only "using K-means". Multiple starts are a setting a
# sensitivity study could examine rather than something to change silently.

#' Starting partition for the multivariate K-centres loop
#'
#' Produces the initial class assignment that \code{\link{MKCFC}} iterates from,
#' implementing Step 2 of the manuscript's algorithm.
#'
#' @section Which start to use: \code{"kmeans"} is the default and is
#' \strong{not} always usable. It clusters one global multivariate FPCA score
#' set, and when that produces a degenerate split the loop aborts on the
#' ten-curve minimum. In simulation on heavy-tailed data it aborted on every
#' replicate, and on real cohorts it can abort outright, while random starts do
#' not. \code{"random"} assigns classes uniformly at random.
#' \code{"two-staged"} tries k-means and falls back to random only when k-means
#' returns a class below ten members or fewer than \code{K} classes.
#'
#' For sensitivity analysis, prefer supplying \code{init.labels} to
#' \code{\link{MKCFC}} directly: strategies such as a hierarchical start or an
#' oracle start cannot be expressed through \code{init.method}.
#'
#' @param Ly A list of \code{P} elements, one per variable, each a list of
#'   \code{n} numeric vectors of observed values.
#' @param Lt A list of the same shape holding observation times.
#' @param K Number of classes.
#' @param init.method One of \code{"kmeans"}, \code{"random"} or
#'   \code{"two-staged"}.
#' @param robust Logical, whether the global FPCA used to build the scores is
#'   the robust one.
#' @param optns Options for the marginal step; see \code{\link{MKCFC}}.
#' @param mfpca_scatter Scatter for the second stage, \code{"sign"} or
#'   \code{"moment"}.
#'
#' @return An integer vector of length \code{n} with values in \code{1:K}.
#'
#' @seealso \code{\link{MKCFC}}, which calls this when \code{init.labels} is not
#'   supplied.
#'
#' @export
MKCFC_init <- function(Ly, Lt, K, init.method, robust, optns,
                       mfpca_scatter = c("sign", "moment")) {
  init.method <- match.arg(init.method, c("kmeans", "random", "two-staged"))
  n <- length(Ly[[1]])

  random_labels <- function()
    sample.int(K, n, replace = TRUE, prob = rep(1 / K, K))

  if (init.method == "random") return(random_labels())

  # MRFPCA when robust, MFPCA otherwise, over ALL subjects as one group, then
  # k-means on the multivariate scores.
  scores <- MFPCA_k(1L, Ly, Lt, rep(1L, n), optns, robust, mfpca_scatter)$multi_score
  clust0 <- stats::kmeans(scores, centers = K)$cluster

  # "or randomly assign the class memberships". With the corrected start this is
  # a genuine last resort rather than the usual path: it fires only when k-means
  # still returns a degenerate split, which is the condition that later makes
  # .mrkcfc_require_cluster_sizes() refuse the fit.
  if (init.method == "two-staged" &&
      (min(table(clust0)) < 10L || length(unique(clust0)) < K)) {
    clust0 <- random_labels()
  }
  clust0
}

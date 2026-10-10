#' Multivariate K-centres functional clustering, with an optional robust
#' two-stage estimator
#'
#' Clusters sparsely and irregularly observed multivariate longitudinal
#' trajectories by iterating between class-specific multivariate functional
#' principal component analysis and reassignment to the best-reconstructing
#' class, after Chiou and Li (2007) extended to several variables.
#'
#' @section The two stages, and which arguments govern which:
#' Estimation proceeds in two stages, and \strong{they are controlled
#' separately}. \code{robust} governs the \emph{marginal} step, choosing Kendall
#' FPCA over ordinary FPCA for each variable; \code{mfpca_scatter} governs the
#' \emph{multivariate} step, choosing a spatial-sign scatter over the sample
#' second moment for the eigendecomposition of the stacked scores.
#'
#' The combination matters more than either setting alone. Setting
#' \code{robust = TRUE} with \code{mfpca_scatter = "moment"} makes the first
#' stage robust and leaves the second stage a moment estimator, and in
#' simulation that pairing recovers the eigenfunction subspace \emph{worse than
#' either consistent choice}: the robust marginal step returns non-negligible
#' variance in components the non-robust one drives to nearly zero, and a moment
#' scatter weights those components while a sign scatter does not.
#'
#' @section Convergence: The loop stops when the reclassification proportion
#' against any earlier iteration falls to \code{conv.prop}. At the default
#' \code{conv.prop = 0} this demands an exact repeat of a previous labelling,
#' which under heavy tails a partition oscillating over a few boundary subjects
#' may never produce. A returned \code{conv = FALSE} therefore does not by
#' itself mean the partition is unsettled; inspect \code{iter} to distinguish a
#' run that exhausted \code{maxIter} from one that stopped early.
#'
#' The loop also aborts if any class falls below ten curves. This is not rare:
#' with \code{init.method = "kmeans"} it occurs on heavy-tailed data and on some
#' real cohorts, while random restarts do not exhibit it. If a fit stops with
#' that error, refit from several random starts and keep the one minimising the
#' reconstruction error rather than lowering the threshold.
#'
#' @param Ly A list of \code{P} elements, one per variable, each itself a list
#'   of \code{n} numeric vectors holding one subject's observed values.
#'   \code{Ly[[j]][[i]]} is subject \code{i}'s measurements on variable
#'   \code{j}. Supply either \code{Ly} and \code{Lt}, or the data-frame form.
#' @param Lt A list of the same shape as \code{Ly} holding the observation
#'   times. \strong{\code{Lt} is a list of \code{P} elements, not of \code{n}}:
#'   the schedule is repeated per variable even when it is shared across them,
#'   and \code{Lt[[j]][[i]]} must have the same length as \code{Ly[[j]][[i]]}.
#' @param data A long-format data frame, as an alternative to \code{Ly} and
#'   \code{Lt}. Must be supplied together with \code{id}, \code{ys} and
#'   \code{timevar}.
#' @param id Name of the subject identifier column in \code{data}.
#' @param ys Character vector naming the \code{P} outcome columns in
#'   \code{data}.
#' @param timevar Name of the observation-time column in \code{data}.
#' @param K Number of classes. A positive integer, at most the number of
#'   subjects.
#' @param robust Logical. \code{TRUE} estimates each variable's marginal
#'   functional principal components by Kendall FPCA, which depends on the
#'   observations only through their directions and so remains defined when
#'   second moments do not exist. \code{FALSE} uses ordinary FPCA.
#' @param optns A list of options for the marginal step. \strong{The two
#'   branches take different lists and reject each other's names.} With
#'   \code{robust = TRUE}: \code{interval}, \code{nK}, \code{bw}, \code{bwK},
#'   \code{bwmean}, \code{nRegGrid} and \code{basis}, the last a
#'   \code{\link[fda]{create.bspline.basis}} object over \code{interval}. With
#'   \code{robust = FALSE}: \code{dataType}, \code{methodMuCovEst},
#'   \code{nRegGrid}, \code{methodSelectK} and \code{useBinnedData}, as taken by
#'   \pkg{fdapace}.
#' @param mfpca_scatter Scatter used in the second stage, for the
#'   eigendecomposition of the stacked univariate scores. \code{"sign"} uses the
#'   spatial-sign covariance matrix, whose eigenvectors coincide with the
#'   covariance's under elliptical symmetry while requiring no moments;
#'   \code{"moment"} uses the sample second moment.
#' @param init.method Starting partition when \code{init.labels} is not given:
#'   \code{"kmeans"} on global multivariate scores, \code{"random"}, or
#'   \code{"two-staged"}. See the note on aborting above before relying on the
#'   default.
#' @param init.labels An integer vector of length \code{n} with values in
#'   \code{1:K}, supplying the starting partition directly and bypassing
#'   \code{init.method}. Intended for studying sensitivity to initialisation,
#'   where strategies such as a hierarchical or oracle start cannot be expressed
#'   through \code{init.method}.
#' @param maxIter Maximum number of reassignment iterations.
#' @param parallel Logical, whether to parallelise the reassignment step. The
#'   parallelism is over \emph{subjects}, so it scales with \code{n} rather than
#'   with \code{K}; the per-class refit remains serial.
#' @param nmcores Number of workers when \code{parallel = TRUE}. Required in
#'   that case, and capped internally at \code{n}.
#' @param trace \code{1} prints per-iteration progress, \code{0} is silent.
#' @param conv.prop Reclassification proportion at or below which the loop is
#'   declared converged. See the convergence section.
#' @param LOO Logical, whether scores are computed leaving each subject out of
#'   its class's estimated components.
#' @param basis_common Logical. \code{FALSE}, the default, estimates
#'   eigenfunctions separately within each class. \code{TRUE} shares one set
#'   across classes and exists for ablation: where the classes genuinely differ
#'   in shape, a common basis is misspecified by construction.
#'
#' @return A list with components
#'   \item{K}{the number of classes fitted.}
#'   \item{res_clust}{a data frame of class assignments, with column \code{g},
#'     and a column \code{id} when the data-frame input form was used.}
#'   \item{init.method}{the starting method recorded as given.}
#'   \item{iter}{iterations performed. \code{0} when \code{K = 1}.}
#'   \item{conv}{logical, whether the reclassification criterion was met.}
#'   \item{optns}{the options list as supplied.}
#'   \item{MFPCA_C}{a list of \code{K} class-specific multivariate FPCA fits, as
#'     returned by \code{\link{MFPCA_k}}, each carrying the eigenfunctions,
#'     multivariate scores and the rotation relating them to the stacked
#'     univariate scores.}
#'
#' @references
#' Chiou, J.-M. and Li, P.-L. (2007). Functional clustering and identifying
#' substructures of longitudinal data. \emph{Journal of the Royal Statistical
#' Society, Series B}, 69, 679--699.
#'
#' Happ, C. and Greven, S. (2018). Multivariate functional principal component
#' analysis for data observed on different (dimensional) domains. \emph{Journal
#' of the American Statistical Association}, 113, 649--659.
#'
#' @seealso \code{\link{MFPCA_k}} for the per-class multivariate step,
#'   \code{\link{MKCFC_init}} for the starting partition, and
#'   \code{\link{silhouette_score}} for choosing \code{K}.
#'
#' @examples
#' \donttest{
#' data(pbc2)
#'
#' # Two markers, subjects with at least two visits on both.
#' vars <- c("serBilir", "albumin")
#' d <- pbc2[stats::complete.cases(pbc2[, c("id", "year", vars)]), ]
#' keep <- names(which(table(d$id) >= 2))
#' d <- d[d$id %in% keep, ]
#'
#' # Time is rescaled to [0, 1] so that the basis and bandwidth below are on the
#' # same footing as the observation grid.
#' rng <- range(d$year)
#' d$t <- (d$year - rng[1]) / diff(rng)
#'
#' ids <- unique(d$id)
#' Ly <- lapply(vars, function(v) lapply(ids, function(i) d[[v]][d$id == i]))
#' Lt <- rep(list(lapply(ids, function(i) d$t[d$id == i])), length(vars))
#'
#' optns <- list(
#'   interval = c(0, 1), nK = 2, bw = 0.1, bwK = "GCV", bwmean = "GCV",
#'   nRegGrid = 51,
#'   basis = fda::create.bspline.basis(c(0, 1), nbasis = 13, norder = 4))
#'
#' fit <- MKCFC(Ly = Ly, Lt = Lt, K = 2, robust = TRUE, optns = optns,
#'              mfpca_scatter = "sign", maxIter = 20, parallel = FALSE,
#'              trace = 0)
#' table(fit$res_clust$g)
#' fit$conv
#' }
#'
#' @export
MKCFC <- function(Ly = NULL, Lt = NULL, data = NULL, id = NULL, ys = NULL,
                  timevar = NULL, K, robust, optns,
                  mfpca_scatter = c("sign", "moment"),
                  init.method = "kmeans", init.labels = NULL, maxIter = 20,
                  parallel = TRUE, nmcores = NULL, trace = 1,
                  conv.prop = 0, LOO = TRUE, basis_common = FALSE) {
  mfpca_scatter <- match.arg(mfpca_scatter)
  list_input <- !is.null(Ly) || !is.null(Lt)
  data_input <- !is.null(data) || !is.null(id) || !is.null(ys) || !is.null(timevar)
  if (list_input == data_input) {
    stop("Supply either Ly and Lt, or data with id, ys, and timevar.", call. = FALSE)
  }
  if (length(K) != 1L || K < 1L || K != as.integer(K)) {
    stop("K must be a positive integer.", call. = FALSE)
  }
  if (length(maxIter) != 1L || maxIter < 1L || maxIter != as.integer(maxIter)) {
    stop("maxIter must be a positive integer.", call. = FALSE)
  }

  uid <- NULL
  if (list_input) {
    if (is.null(Ly) || is.null(Lt)) {
      stop("Both Ly and Lt must be supplied together.", call. = FALSE)
    }
    n <- check_input(Ly, Lt)
  } else {
    if (is.null(data) || is.null(id) || is.null(ys) || is.null(timevar)) {
      stop("data, id, ys, and timevar must be supplied together.", call. = FALSE)
    }
    dt <- MFPCAInputs(id, ys, timevar, data)
    Ly <- dt$Ly
    Lt <- dt$Lt
    uid <- unlist(dt$id)
    n <- length(uid)
  }
  if (K > n) {
    stop("K cannot exceed the number of subjects.", call. = FALSE)
  }

  workers <- .mrkcfc_validate_parallel_workers(parallel, nmcores, n)
  cluster <- NULL
  if (workers > 1L) {
    cluster <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    .mrkcfc_prepare_workers(
      cluster,
      c("clust_update", "predicted_yhat", ".mrkcfc_solve_covariance",
        ".mrkcfc_nonrobust_scores", ".mrkcfc_mfpca_rotation",
        ".mrkcfc_mfpca_reconstruct", "robust_stage2", "sscm", "spatial_median")
    )
  }

  if (K == 1L) {
    final_labels <- .mrkcfc_one_cluster_labels(n)
    MFPCA_C <- lapply(seq_len(K), function(k) {
      MFPCA_k(k, Ly, Lt, final_labels, optns, robust, mfpca_scatter)
    })
    res_clust <- if (is.null(uid)) {
      data.frame(g = final_labels)
    } else {
      data.frame(id = uid, g = final_labels)
    }
    return(list(
      K = K, res_clust = res_clust, init.method = init.method,
      iter = 0L, conv = TRUE, optns = optns, MFPCA_C = MFPCA_C
    ))
  }

  # init.labels supplies a starting partition directly, bypassing MKCFC_init().
  # An initialisation study needs it: comparing the k-means start against a
  # hierarchical or an oracle true-label start is not expressible through
  # init.method. It is also the honest way to study initialisation sensitivity,
  # since the alternative is reimplementing the loop outside the package and
  # then arguing that the reimplementation matches.
  clust0 <- if (is.null(init.labels)) {
    MKCFC_init(Ly, Lt, K, init.method, robust, optns, mfpca_scatter)
  } else {
    lab <- as.integer(init.labels)
    if (length(lab) != n || anyNA(lab) || !all(lab %in% seq_len(K)))
      stop("init.labels must be a vector of length n with values in 1..K.",
           call. = FALSE)
    lab
  }
  .mrkcfc_require_cluster_sizes(clust0, K)
  clustlist <- list(clust0)
  final_labels <- clust0
  iteration <- 1L
  conv <- FALSE

  if (maxIter >= 2L) {
    for (candidate_iteration in seq.int(2L, as.integer(maxIter))) {
      if (trace == 2) {
        cat(candidate_iteration, "th iteration\n")
      }
      ClustIds <- .mrkcfc_cluster_ids(final_labels, K)
      # basis_common = TRUE is the class-common ablation: one set of eigenfunctions
      # for every cluster, per-cluster means only. It changes the REFERENCE the
      # reassignment step scores against, which is what makes it a change of
      # method rather than of reporting.
      ref_hatyc <- .mrkcfc_fit_multivariate_reference(
        Ly, Lt, K, final_labels, optns, robust, basis_common
      )
      if (is.null(cluster)) {
        clust_it <- lapply(seq_len(n), clust_update,
                           ClustIds = ClustIds, Ly = Ly, Lt = Lt, K = K,
                           ref_hatyc = ref_hatyc, optns = optns,
                           robust = robust, LOO = LOO,
                           mfpca_scatter = mfpca_scatter)
      } else {
        clust_it <- parallel::parLapply(
          cluster, seq_len(n), clust_update,
          ClustIds = ClustIds, Ly = Ly, Lt = Lt, K = K,
          ref_hatyc = ref_hatyc, optns = optns, robust = robust, LOO = LOO,
          mfpca_scatter = mfpca_scatter
        )
      }
      candidate <- vapply(clust_it, function(x) x$clust_i, integer(1))
      candidate_ids <- .mrkcfc_cluster_ids(candidate, K)
      if (any(lengths(candidate_ids) < 10L)) {
        if (trace >= 1) cat("Error: Less than 10 curves in a group.\n")
        break
      }
      clustlist[[candidate_iteration]] <- candidate
      final_labels <- candidate
      iteration <- candidate_iteration
      converged <- any(vapply(clustlist[-length(clustlist)], function(previous) {
        mclust::classError(previous, candidate)$errorRate <= conv.prop
      }, logical(1)))
      if (converged) {
        if (trace >= 1) cat("Clustering algorithm converged.\n")
        conv <- TRUE
        break
      }
    }
  }

  MFPCA_C <- lapply(seq_len(K), function(k) {
    MFPCA_k(k, Ly, Lt, final_labels, optns, robust, mfpca_scatter)
  })
  res_clust <- if (is.null(uid)) {
    data.frame(g = final_labels)
  } else {
    data.frame(id = uid, g = final_labels)
  }
  list(
    K = K, res_clust = res_clust, init.method = init.method,
    iter = iteration, conv = conv, optns = optns, MFPCA_C = MFPCA_C
  )
}

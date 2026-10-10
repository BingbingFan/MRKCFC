.mrkcfc_validate_univariate <- function(Ly, Lt) {
  if (!is.list(Ly) || !is.list(Lt) || length(Ly) == 0L || length(Lt) == 0L) {
    stop("Ly and Lt must be non-empty lists.", call. = FALSE)
  }
  if (length(Ly) != length(Lt)) {
    stop("Ly and Lt must contain the same number of subjects.", call. = FALSE)
  }
  valid <- vapply(seq_along(Ly), function(i) {
    length(Ly[[i]]) > 0L && length(Ly[[i]]) == length(Lt[[i]])
  }, logical(1))
  if (!all(valid)) {
    stop("Each subject's Ly and Lt vectors must be non-empty and have the same length.", call. = FALSE)
  }
  length(Ly)
}

.mrkcfc_validate_multivariate <- function(Ly, Lt) {
  if (!is.list(Ly) || !is.list(Lt) || length(Ly) == 0L || length(Lt) == 0L) {
    stop("Ly and Lt must be non-empty lists of outcomes.", call. = FALSE)
  }
  if (length(Ly) != length(Lt)) {
    stop("Ly and Lt must contain the same number of outcomes.", call. = FALSE)
  }
  n <- .mrkcfc_validate_univariate(Ly[[1]], Lt[[1]])
  n_by_outcome <- vapply(seq_along(Ly), function(j) {
    .mrkcfc_validate_univariate(Ly[[j]], Lt[[j]])
  }, integer(1))
  if (!all(n_by_outcome == n)) {
    stop("Every outcome must contain the same number of subjects.", call. = FALSE)
  }
  n
}

.mrkcfc_cluster_ids <- function(labels, K) {
  if (length(labels) == 0L || length(K) != 1L || K < 1L || K != as.integer(K)) {
    stop("K must be a positive integer and labels must be non-empty.", call. = FALSE)
  }
  lapply(seq_len(as.integer(K)), function(k) which(labels == k))
}

.mrkcfc_require_cluster_sizes <- function(labels, K, minimum = 10L) {
  ids <- .mrkcfc_cluster_ids(labels, K)
  if (any(lengths(ids) < minimum)) {
    stop(sprintf("Each cluster must contain at least %d curves.", minimum), call. = FALSE)
  }
  invisible(ids)
}

.mrkcfc_one_cluster_labels <- function(n) {
  if (length(n) != 1L || n < 1L || n != as.integer(n)) {
    stop("n must be a positive integer.", call. = FALSE)
  }
  rep.int(1L, as.integer(n))
}

.mrkcfc_validate_parallel_workers <- function(parallel, nmcores, n) {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("parallel must be TRUE or FALSE.", call. = FALSE)
  }
  if (!parallel) {
    return(1L)
  }
  if (is.null(nmcores) || length(nmcores) != 1L || !is.finite(nmcores) ||
      nmcores < 1L || nmcores != as.integer(nmcores)) {
    stop("nmcores must be a positive integer when parallel = TRUE.", call. = FALSE)
  }
  min(as.integer(nmcores), as.integer(n))
}

# Make the workers able to run the update functions. Two modes must both work,
# and they fail differently, which is why this is not a plain clusterExport().
#
#   installed  The update functions live in the MRKCFC namespace. When a worker
#              deserialises one, it resolves that namespace from the WORKER's
#              own .libPaths(). If a different build of MRKCFC is installed
#              there, the worker silently runs THAT one. The failure surfaces
#              as "unused argument", not "object not found", because the worker
#              found a real function of the same name from an older version.
#              Pushing the master's library path fixes it, and loading the
#              namespace explicitly makes the dependency visible rather than
#              relying on deserialisation to trigger it.
#
#   sourced    R/*.R has been source()d, so the functions are in the global
#              environment and there may be no installed MRKCFC at all, or only
#              a stale one. Here the closures themselves must be shipped, and
#              loading a namespace would be actively wrong.
#
# The old code always took the second branch, so it worked when sourced and
# broke once installed. Keep both paths; a hand-maintained export list is only
# correct for the sourced case, and forgetting a name there fails solely under
# parallel = TRUE.
.mrkcfc_prepare_workers <- function(cluster, exports) {
  env <- topenv(environment(.mrkcfc_prepare_workers))
  if (isNamespace(env)) {
    parallel::clusterCall(cluster, function(paths) .libPaths(paths), .libPaths())
    pkg <- environmentName(env)
    parallel::clusterCall(cluster, function(p) {
      suppressMessages(requireNamespace(p, quietly = TRUE))
      invisible(NULL)
    }, pkg)
    return(invisible("namespace"))
  }
  parallel::clusterExport(cluster, exports, envir = env)
  invisible("global")
}

.mrkcfc_fit_multivariate_reference <- function(Ly, Lt, K, labels, optns, robust,
                                              basis_common = FALSE) {
  cluster_ids <- .mrkcfc_cluster_ids(labels, K)
  fit_one <- function(ids, j) {
    yy <- lapply(ids, function(i) Ly[[j]][[i]])
    tt <- lapply(ids, function(i) Lt[[j]][[i]])
    if (robust) {
      KFPCA(
        Lt = tt, Ly = yy, interval = optns$interval, nK = optns$nK,
        bw = optns$bw, bwK = optns$bwK, nRegGrid = optns$nRegGrid,
        fdParobj = optns$basis
      )
    } else {
      fdapace::FPCA(yy, tt, optns = optns)
    }
  }

  if (!isTRUE(basis_common)) {
    return(lapply(seq_len(K), function(k)
      lapply(seq_along(Ly), function(j) fit_one(cluster_ids[[k]], j))))
  }

  # --- class-common basis, the ablation ---------------------------------------
  #
  # One set of eigenfunctions estimated from ALL subjects, shared by every
  # cluster; only the MEAN stays per-cluster. If the mean were shared too the
  # clusters would be indistinguishable, every seCost would be equal, and the
  # reassignment step would pick arbitrarily.
  #
  # This is the arm that measures what class-specific basis estimation buys. It
  # is misspecified by construction on this design, where variables 1 and 2 have
  # sines for class 1 and cosines for class 2 -- which is the point: the size of
  # the loss is the result.
  #
  # WHAT HAS TO BE SWAPPED, AND WHY IT IS NOT JUST $mean. predict.KFPCA() does
  # NOT read object$mean. It re-estimates the mean from object$Lt and object$Ly,
  # the fit's own training data (see predict.KFPCA.R, `meanest_new <-
  # MeanEst(Lt, Ly, ...)`). Swapping only $mean would therefore leave the robust
  # branch silently centring on the global mean while appearing to be
  # per-cluster. The training data must be swapped as well.
  #
  # The non-robust branch has no such hidden path: predicted_yhat() forms the
  # residual as `ith_yy - mu_ij` from $mu directly, so swapping $mu suffices.
  #
  # $score and $xiEst are deliberately left at their global values. They are
  # used only for n_M and for the stage-two rotation, and a rotation common
  # across clusters is the correct object for a common-basis arm.
  global <- lapply(seq_along(Ly), function(j) fit_one(seq_along(Ly[[j]]), j))

  lapply(seq_len(K), function(k) {
    ids <- cluster_ids[[k]]
    lapply(seq_along(Ly), function(j) {
      ref <- global[[j]]
      loc <- fit_one(ids, j)          # fitted for its MEAN only
      if (robust) {
        ref$mean <- loc$mean
        ref$Lt <- lapply(ids, function(i) Lt[[j]][[i]])
        ref$Ly <- lapply(ids, function(i) Ly[[j]][[i]])
      } else {
        ref$mu <- loc$mu
      }
      ref
    })
  })
}

.mrkcfc_solve_covariance <- function(covariance, residual,
                                    tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(covariance) || nrow(covariance) != ncol(covariance)) {
    stop("covariance must be a square matrix.", call. = FALSE)
  }
  if (length(residual) != nrow(covariance)) {
    stop("residual must have one value per covariance row.", call. = FALSE)
  }
  covariance <- (covariance + t(covariance)) / 2
  direct <- tryCatch(solve(covariance, residual), error = function(e) NULL)
  if (!is.null(direct)) {
    return(direct)
  }
  scale <- max(1, max(abs(diag(covariance))))
  ridge <- sqrt(.Machine$double.eps) * scale
  regularized <- tryCatch(
    solve(covariance + diag(ridge, nrow(covariance)), residual),
    error = function(e) NULL
  )
  if (!is.null(regularized)) {
    return(regularized)
  }
  MASS::ginv(covariance, tol = tol) %*% residual
}

.mrkcfc_nonrobust_scores <- function(phi, lambda, sigma2, residual) {
  lambda <- pmax(as.numeric(lambda), 0)
  if (!is.matrix(phi) || ncol(phi) != length(lambda)) {
    stop("phi columns must match the number of eigenvalues.", call. = FALSE)
  }
  if (length(residual) != nrow(phi)) {
    stop("residual length must match the number of observation times.", call. = FALSE)
  }
  noise <- max(as.numeric(sigma2), .Machine$double.eps)
  weighted_phi <- sweep(phi, 2L, sqrt(lambda), "*")
  covariance <- tcrossprod(weighted_phi) + diag(noise, nrow(phi))
  as.vector(lambda * crossprod(phi, .mrkcfc_solve_covariance(covariance, residual)))
}

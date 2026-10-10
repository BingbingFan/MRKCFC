ch_index         <- MRKCFC:::ch_index
silhouette_index <- MRKCFC:::silhouette_index
whiten_scores    <- MRKCFC:::whiten_scores

# Two well-separated Gaussian clusters in four dimensions, and nested
# partitions that subdivide the truth. Subdividing a correct cluster must be
# penalised by a usable criterion.
set.seed(11)
n <- 200
truth <- rep(1:2, each = n / 2)
X <- matrix(rnorm(n * 4), n, 4) + cbind(truth * 4, 0, 0, 0)
parts <- list(truth,
              truth + 2L * (seq_len(n) %% 2L),
              truth + 2L * (seq_len(n) %% 3L))

# the bare between/within ratio, without the degrees-of-freedom divisors
raw_ratio <- function(mat, g) {
  g <- as.integer(as.factor(g)); m <- colMeans(mat); B <- 0; W <- 0
  for (k in unique(g)) {
    Xk <- mat[g == k, , drop = FALSE]; mk <- colMeans(Xk)
    B <- B + nrow(Xk) * sum((mk - m)^2)
    W <- W + sum((Xk - rep(mk, each = nrow(Xk)))^2)
  }
  B / W
}

test_that("the Calinski-Harabasz index is maximised at the true K", {
  ch  <- vapply(parts, ch_index, numeric(1), scores = X, robust = FALSE)
  raw <- vapply(parts, raw_ratio, numeric(1), mat = X)
  expect_equal(which.max(ch), 1L)
  # the bare ratio increases with K, so it could only ever select the largest K
  expect_true(all(diff(raw) > 0))
})

test_that("the silhouette is maximised at the true K and stays in [-1, 1]", {
  sil <- vapply(parts, silhouette_index, numeric(1), scores = X, robust = FALSE)
  expect_equal(which.max(sil), 1L)
  expect_true(all(sil >= -1 & sil <= 1))
})

test_that("the silhouette follows Rousseeuw's definition, computed by hand", {
  # a(i) excludes i itself: for a member of a three-point cluster it is the mean
  # of the two distances to its neighbours, not of three distances one of which
  # is zero. A singleton cluster has s(i) = 0 by definition.
  Y <- rbind(c(0, 0), c(1, 0), c(0, 1), c(50, 50))
  g <- c(1L, 1L, 1L, 2L)
  D <- as.matrix(dist(Y))
  s <- function(i) {
    a <- mean(D[i, setdiff(1:3, i)]); b <- D[i, 4]
    (b - a) / max(a, b)
  }
  manual <- mean(c(s(1), s(2), s(3), 0))
  expect_equal(silhouette_index(Y, g, robust = FALSE), manual)
})

test_that("the Mahalanobis whitening removes the effect of an irrelevant scale", {
  Z <- X; Z[, 4] <- Z[, 4] * 100
  ch_m <- ch_index(X, truth, "mahalanobis", robust = FALSE)
  expect_lt(abs(ch_index(Z, truth, "mahalanobis", robust = FALSE) - ch_m), 0.05 * ch_m)
  expect_lt(ch_index(Z, truth, "euclidean", robust = FALSE),
            0.5 * ch_index(X, truth, "euclidean", robust = FALSE))
  W4 <- whiten_scores(X, robust = FALSE)
  expect_lt(max(abs(cov(W4) - diag(4))), 0.15)
})

test_that("the robust criteria are finite under infinite variance", {
  set.seed(4)
  Xt <- matrix(rt(n * 4, df = 1.5), n, 4) + cbind(truth * 6, 0, 0, 0)
  ch_r  <- ch_index(Xt, truth, "mahalanobis", robust = TRUE)
  sil_r <- silhouette_index(Xt, truth, "mahalanobis", robust = TRUE)
  expect_true(is.finite(ch_r) && ch_r > 0)
  expect_true(is.finite(sil_r) && sil_r > 0)
})

test_that("whitening by the within-cluster scatter would make W constant in K", {
  # which is why whiten_scores() uses a global scatter
  wss_within_metric <- function(mat, g) {
    g <- as.integer(as.factor(g))
    resid <- do.call(rbind, lapply(unique(g), function(k) {
      Xk <- mat[g == k, , drop = FALSE]
      Xk - rep(colMeans(Xk), each = nrow(Xk))
    }))
    S <- crossprod(resid) / nrow(resid)
    sum(diag(resid %*% solve(S) %*% t(resid)))
  }
  w_by_K <- vapply(parts, wss_within_metric, numeric(1), mat = X)
  expect_lt(max(abs(w_by_K - n * 4)), 1e-6)
})

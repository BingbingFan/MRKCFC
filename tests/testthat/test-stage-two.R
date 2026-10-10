rotation_of   <- MRKCFC:::.mrkcfc_mfpca_rotation
reconstruct   <- MRKCFC:::.mrkcfc_mfpca_reconstruct
robust_stage2 <- MRKCFC:::robust_stage2
gen_scores    <- MRKCFC:::gen_scores

set.seed(826)
X <- matrix(rnorm(200 * 9), 200, 9) %*% diag(c(5, 4, 3.5, 3, 2.5, 2, 1.5, 1, 0.5))

test_that("the moment arm is the symmetric eigendecomposition of the second-moment matrix", {
  expect_identical(rotation_of(X, "moment"), robust_stage2(X, robust = FALSE)$vectors)
  old_svd <- svd(crossprod(X) / (nrow(X) - 1))$u
  # the svd of the same matrix spans the same subspace; the columns differ only by sign
  expect_equal(abs(rotation_of(X, "moment")), abs(old_svd))
})

test_that("the sign arm is the spatial-sign rotation, and it is the default", {
  expect_identical(rotation_of(X, "sign"), robust_stage2(X, robust = TRUE)$vectors)
  old_svd <- svd(crossprod(X) / (nrow(X) - 1))$u
  expect_gt(max(abs(abs(rotation_of(X, "sign")) - abs(old_svd))), 1e-6)
  expect_identical(rotation_of(X), rotation_of(X, "sign"))
  expect_error(rotation_of(X, "banana"))
})

test_that("the reconstruction uses one score per component from the full stacked vector", {
  set.seed(99)
  nM   <- list(3L, 3L, 3L)
  nobs <- c(5L, 4L, 6L)
  pi_obs <- lapply(seq_len(3), function(j) matrix(rnorm(nobs[j] * nM[[j]]), nobs[j], nM[[j]]))
  xi_i   <- lapply(nM, function(m) rnorm(m))
  mu_obs <- lapply(nobs, function(n) rnorm(n))
  rot    <- rotation_of(X, "moment")

  r0 <- reconstruct(rot, pi_obs, xi_i, nM, mu_obs)
  expect_equal(r0$multi_score, as.vector(crossprod(rot, unlist(xi_i, use.names = FALSE))))
  expect_length(r0$multi_score, nrow(rot))
  expect_equal(t(r0$multi_eigen[[1]]), pi_obs[[1]] %*% rot[1:3, 1:3, drop = FALSE])

  # perturbing variable 3's scores must move variable 1's reconstruction
  xi_pert <- xi_i; xi_pert[[3]] <- xi_pert[[3]] + 1
  r1 <- reconstruct(rot, pi_obs, xi_pert, nM, mu_obs)
  expect_gt(max(abs(r1$fittedX[[1]] - r0$fittedX[[1]])), 1e-8)

  # invariant to the sign of the eigenvectors, and so to the solver that produced them
  D  <- diag(sample(c(-1, 1), nrow(rot), TRUE))
  rD <- reconstruct(rot %*% D, pi_obs, xi_i, nM, mu_obs)
  expect_equal(rD$fittedX, r0$fittedX)
  old_svd <- svd(crossprod(X) / (nrow(X) - 1))$u
  expect_equal(reconstruct(old_svd, pi_obs, xi_i, nM, mu_obs)$fittedX, r0$fittedX)

  # the truncation binds: retaining every component would collapse to the univariate fit
  univariate <- as.numeric(mu_obs[[1]] + pi_obs[[1]] %*% xi_i[[1]])
  expect_gt(max(abs(r0$fittedX[[1]] - univariate)), 1e-8)
})

test_that("mfpca_scatter is an argument of every function on the fitting path", {
  for (fn in c("MKCFC", "MFPCA_k", "clust_update", "predicted_yhat"))
    expect_true("mfpca_scatter" %in% names(formals(get(fn, envir = asNamespace("MRKCFC")))))
})

test_that("the sign arm orders components by the recovered robust spread", {
  set.seed(827)
  Zc <- matrix(0.2, 6, 6); bl <- list(c(9, 4), c(8, 3.5), c(8.5, 3.5))
  for (b in 0:2) Zc[(2 * b + 1):(2 * b + 2), (2 * b + 1):(2 * b + 2)] <- diag(bl[[b + 1]])
  for (arm in c(1L, 4L, 6L)) {
    Xa <- gen_scores(2000, arm, Zc)
    for (rob in c(FALSE, TRUE)) {
      v <- robust_stage2(Xa, robust = rob)$values
      expect_false(is.unsorted(rev(v)))
    }
    s <- robust_stage2(Xa, robust = TRUE)
    spread <- apply(Xa %*% s$vectors, 2, function(z) (IQR(z) / 1.349)^2)
    expect_equal(spread, s$values)
  }
  X6 <- gen_scores(2000, 6L, Zc)
  expect_equal(robust_stage2(X6, M_keep = 3L)$values, head(robust_stage2(X6)$values, 3))
})

test_that("the starting partition reads every variable", {
  expect_true(all(c("robust", "mfpca_scatter") %in% names(formals(MKCFC_init))))
  set.seed(4242)
  ni <- lapply(seq_len(40), function(i) sort(runif(sample(4:9, 1))))
  mk <- function(shift3 = 0) {
    Ly <- lapply(1:3, function(p) lapply(seq_along(ni), function(i)
      sin(2 * pi * ni[[i]]) * p + rnorm(length(ni[[i]]), 0, 0.2) +
        (if (p == 3L) shift3 * (i %% 2) else 0)))
    list(Ly = Ly, Lt = lapply(1:3, function(p) ni))
  }
  op <- list(dataType = "Sparse", methodMuCovEst = "smooth", nRegGrid = 21,
             methodSelectK = 2, useBinnedData = "OFF")
  d0 <- mk(0); d3 <- mk(6)
  set.seed(1); i0 <- MKCFC_init(d0$Ly, d0$Lt, 2L, "kmeans", FALSE, op, "moment")
  set.seed(1); i3 <- MKCFC_init(d3$Ly, d3$Lt, 2L, "kmeans", FALSE, op, "moment")
  expect_false(identical(as.integer(i0), as.integer(i3)))
  expect_error(MKCFC_init(d0$Ly, d0$Lt, 2L, "banana", FALSE, op, "moment"))
  expect_length(MKCFC_init(d0$Ly, d0$Lt, 2L, "random", FALSE, op, "moment"), length(ni))
})

test_that("parallel workers are prepared through the helper", {
  prep <- MRKCFC:::.mrkcfc_prepare_workers
  expect_true(is.function(prep))
  expect_true(all(c("isNamespace", "clusterCall", "clusterExport") %in% all.names(body(prep))))
  expect_true(".mrkcfc_prepare_workers" %in% all.names(body(MKCFC)))
})

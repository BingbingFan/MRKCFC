# The class-common basis option: every class shares one set of eigenfunctions
# while the means stay per class. On the robust branch the training data must
# be swapped as well, because predict.KFPCA() re-estimates the mean from the
# object's own Lt and Ly rather than reading its mean component.

fit_ref <- MRKCFC:::.mrkcfc_fit_multivariate_reference

make_data <- function() {
  set.seed(7)
  n <- 60; J <- 2; tt <- seq(0, 1, length.out = 101)
  labels <- rep(1:2, each = n / 2)
  Ly <- Lt <- vector("list", J)
  for (j in seq_len(J)) {
    Ly[[j]] <- Lt[[j]] <- vector("list", n)
    for (i in seq_len(n)) {
      ti <- sort(sample(tt[-c(1, length(tt))], 7))
      # class 1 sines, class 2 cosines: a class-specific shape, so a common
      # basis is misspecified here
      f <- if (labels[i] == 1) 2 * sin(2 * pi * ti) else 2 * cos(2 * pi * ti)
      Ly[[j]][[i]] <- f + 3 * (labels[i] == 1) + rnorm(7, 0, .3)
      Lt[[j]][[i]] <- ti
    }
  }
  list(Ly = Ly, Lt = Lt, labels = labels, n = n)
}

test_that("a common basis shares eigenfunctions and keeps per-class means (non-robust)", {
  skip_on_cran()
  d <- make_data()
  optns <- list(dataType = "Sparse", methodMuCovEst = "smooth",
                nRegGrid = 51, methodSelectK = 2L, useBinnedData = "OFF")
  spec <- fit_ref(d$Ly, d$Lt, 2L, d$labels, optns, robust = FALSE, basis_common = FALSE)
  comm <- fit_ref(d$Ly, d$Lt, 2L, d$labels, optns, robust = FALSE, basis_common = TRUE)
  expect_false(isTRUE(all.equal(spec[[1]][[1]]$phi, spec[[2]][[1]]$phi)))
  expect_equal(comm[[1]][[1]]$phi, comm[[2]][[1]]$phi)
  expect_false(isTRUE(all.equal(comm[[1]][[1]]$mu, comm[[2]][[1]]$mu)))
  # the class mean is the class's own fitted mean, not the global one
  loc1 <- fit_ref(d$Ly, d$Lt, 2L, d$labels, optns, robust = FALSE)[[1]][[1]]
  expect_equal(comm[[1]][[1]]$mu, loc1$mu)
})

test_that("on the robust branch the training data are swapped with the mean", {
  skip_on_cran()
  d <- make_data()
  roptns <- list(interval = c(0, 1), nK = 2L, bw = 0.1, bwK = "GCV", nRegGrid = 51,
                 basis = fda::create.bspline.basis(c(0, 1), nbasis = 8))
  rcomm <- fit_ref(d$Ly, d$Lt, 2L, d$labels, roptns, robust = TRUE, basis_common = TRUE)
  expect_equal(rcomm[[1]][[1]]$FPC_dis, rcomm[[2]][[1]]$FPC_dis)
  expect_false(isTRUE(all.equal(rcomm[[1]][[1]]$mean, rcomm[[2]][[1]]$mean)))
  expect_equal(length(rcomm[[1]][[1]]$Lt), d$n / 2)
  expect_equal(length(rcomm[[2]][[1]]$Lt), d$n / 2)
  expect_false(identical(rcomm[[1]][[1]]$Ly, rcomm[[2]][[1]]$Ly))
  # and predict() centres differently per class
  p1 <- stats::predict(rcomm[[1]][[1]], newLt = list(d$Lt[[1]][[1]]),
                       newLy = list(d$Ly[[1]][[1]]), nK = 2L, more = TRUE)
  p2 <- stats::predict(rcomm[[2]][[1]], newLt = list(d$Lt[[1]][[1]]),
                       newLy = list(d$Ly[[1]][[1]]), nK = 2L, more = TRUE)
  expect_false(isTRUE(all.equal(p1$score_new, p2$score_new)))
})

test_that("MKCFC() returns a partition with the common basis", {
  skip_on_cran()
  d <- make_data()
  optns <- list(dataType = "Sparse", methodMuCovEst = "smooth",
                nRegGrid = 51, methodSelectK = 2L, useBinnedData = "OFF")
  f1 <- MKCFC(Ly = d$Ly, Lt = d$Lt, K = 2L, robust = FALSE, optns = optns,
              mfpca_scatter = "moment", init.method = "two-staged",
              maxIter = 5L, parallel = FALSE, trace = 0, LOO = FALSE,
              basis_common = TRUE)
  g1 <- f1$res_clust$g
  expect_length(g1, d$n)
  expect_true(all(g1 %in% 1:2))
  expect_gte(min(table(g1)), 1L)
})

test_that("the covariance solver returns a direct solution", {
  result <- MRKCFC:::.mrkcfc_solve_covariance(diag(2), c(3, -1))
  expect_equal(as.vector(result), c(3, -1))
})

test_that("the covariance solver handles a singular matrix", {
  result <- MRKCFC:::.mrkcfc_solve_covariance(matrix(1, 2, 2), c(2, 2))
  expect_equal(as.vector(result), c(1, 1), tolerance = 1e-7)
})

test_that("non-robust scores use the model covariance", {
  score <- MRKCFC:::.mrkcfc_nonrobust_scores(diag(2), c(4, 1), 1, c(2, 2))
  expect_equal(score, c(8 / 5, 1))
})

test_that("the class index and size helpers behave", {
  expect_equal(MRKCFC:::.mrkcfc_cluster_ids(c(1L, 2L, 1L), 2L), list(c(1L, 3L), 2L))
  expect_equal(MRKCFC:::.mrkcfc_one_cluster_labels(4L), rep(1L, 4L))
  expect_error(MRKCFC:::.mrkcfc_require_cluster_sizes(c(1L, 1L, 2L), 2L), "at least 10")
  expect_error(MRKCFC:::.mrkcfc_validate_univariate(list(1:2), list(1)), "same length")
})

test_that("the long-format input is split by subject in order of appearance", {
  input_data <- data.frame(subject = c("b", "b", "a"), time = c(2, 1, 1), y = c(20, 10, 30))
  input_result <- MFPCAInputs("subject", "y", "time", input_data)
  expect_equal(unlist(input_result$id), c("b", "a"))
  expect_equal(input_result$Lt[[1]][[1]], c(1, 2))
  expect_equal(input_result$Ly[[1]][[1]], c(10, 20))
})

test_that("the KFPCA runtime is internal to the package", {
  # This reads the source tree, which exists in a checkout but not under
  # R CMD check, where the tests run against an installed copy; skip there.
  package_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
                                mustWork = FALSE)
  testthat::skip_if_not(
    file.exists(file.path(package_root, "DESCRIPTION")) &&
      dir.exists(file.path(package_root, "R")),
    "source tree not available; this check applies to the sources, not an install")
  description_lines <- readLines(file.path(package_root, "DESCRIPTION"), warn = FALSE)
  source_paths <- list.files(file.path(package_root, "R"), pattern = "\\.R$", full.names = TRUE)
  source_text <- unlist(lapply(source_paths, readLines, warn = FALSE))

  expect_false(any(grepl("^\\s*KFPCA\\s*(,|$)", description_lines)))
  expect_false(any(grepl("KFPCA::", source_text, fixed = TRUE)))
  expect_true(is.function(MRKCFC::KFPCA))
  expect_true(is.function(MRKCFC:::predict.KFPCA))
  expect_match(
    paste(readLines(file.path(package_root, "NAMESPACE"), warn = FALSE), collapse = "\\n"),
    "S3method\\(predict,KFPCA\\)"
  )
})

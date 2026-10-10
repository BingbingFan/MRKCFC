#' @title Relative reconstruction error
#'
#' @description
#' For each variable, the residual sum of squares divided by the total sum of
#' squares about that variable's mean: 0 is perfect reconstruction, 1 is no
#' better than the mean of the observations, and values above 1 mean the fit
#' reconstructs worse than a constant. The overall value is the mean over
#' variables.
#'
#' @details
#' A sum of squared errors across variables adds quantities in different units,
#' and is coherent only because normalization precedes it; it also depends on
#' the scale the data are reported on. Relativizing per variable makes the
#' measure unit-free, comparable between a simulation and an application, and
#' readable on its own.
#'
#' The quantity is an \emph{in-sample} reconstruction error whenever
#' \code{yhat} comes from a fit that used \code{y}. For predictive performance,
#' withhold observations from the fit and pass only those here; in simulation,
#' pass the true latent trajectory as \code{y}.
#'
#' \code{trim} additionally restricts each variable to its middle
#' \code{1 - trim} of observations, by absolute deviation from the median. If
#' the trimmed and untrimmed values agree, an advantage is broad-based; if they
#' diverge, it is concentrated in the tail, which is a different claim.
#' Trimming is on how extreme the observation is, not on how large the error
#' is: trimming by error would remove each method's own worst cases and score
#' two methods on different subsets.
#'
#' @param y Numeric vector of observed values, all variables stacked.
#' @param yhat Numeric vector of reconstructed values, the same length as
#'   \code{y}.
#' @param variable A vector the same length as \code{y} giving each
#'   observation's variable.
#' @param trim Proportion of the most extreme observations of each variable to
#'   exclude. \code{0} uses all observations.
#'
#' @return A named numeric vector: one \code{per_variable} entry per variable
#'   and \code{overall}, their mean.
#'
#' @seealso \code{\link{predicted_yhat}}, \code{\link{robust_normalize}}
#'
#' @examples
#' set.seed(1)
#' y <- stats::rnorm(100)
#' v <- rep(1:2, each = 50)
#' rre(y, yhat = 0.8 * y, variable = v)              # reconstructs most of the variation
#' rre(y, yhat = ave(y, v), variable = v)            # exactly 1: no better than the mean
#'
#' @export
rre <- function(y, yhat, variable, trim = 0) {
  stopifnot(length(y) == length(yhat), length(y) == length(variable))
  per_var <- tapply(seq_along(y), variable, function(i) {
    if (trim > 0) {
      dev  <- abs(y[i] - median(y[i]))
      i <- i[dev <= quantile(dev, 1 - trim, names = FALSE)]
    }
    sum((y[i] - yhat[i])^2) / sum((y[i] - mean(y[i]))^2)
  })
  c(per_variable = per_var, overall = mean(per_var))
}

# Guard: assay measurements have mass at the detection floor. If more than
# three quarters of the observations of a variable sit there, Q1 == Q3 and the
# IQR is zero. With type-7 quantiles the exact condition is f >= 0.75n + 0.25.
# Half is not enough: between 50 and 75 percent the median sits on the floor
# but Q3 does not, so the IQR stays positive and this guard does not fire.
# The fallback widens the quantile spread rather than failing.

#' @title Rescale a variable by its median and interquartile range
#'
#' @description
#' Centres at the median and scales by the interquartile range, pooled over all
#' observations of a variable across subjects and time points. The result is
#' centred at 0 with unit interquartile range.
#'
#' @details
#' Min-max rescaling, the obvious alternative, makes the scale of every subject
#' depend on a single extreme observation, which is exactly the sensitivity a
#' robust method exists to remove. The median and interquartile range serve the
#' same purpose of placing a variable on a known scale while depending only on
#' the bulk of the data.
#'
#' Pool, do not standardise within subject. Between-subject differences in level
#' are what the clustering uses, and normalising within subject would delete the
#' signal.
#'
#' Centring is harmless for \code{\link{MKCFC}}: a location shift applies to the
#' observations and their reconstruction alike and cancels in the assignment
#' rule, and the mean function is removed by the FPCA step regardless.
#'
#' The interquartile range is preferred to the median absolute deviation here.
#' The MAD symmetrises deviations about the median, which for a strongly
#' right-skewed marker is governed largely by the left side and can understate
#' the spread. A spread between two order statistics of the bulk behaves better.
#'
#' @section Guard against a zero denominator:
#' Assay measurements often have mass at a detection floor. With type-7
#' quantiles, \eqn{Q_1 = Q_3} requires roughly three quarters of the
#' observations to sit at the floor, not half: between 50 and 75 percent the
#' median is on the floor but \eqn{Q_3} is not, so the interquartile range stays
#' positive. When it does reach zero the function falls back to a wider quantile
#' spread and warns, rather than failing.
#'
#' @param y Numeric vector of all observations of one variable.
#' @param fallback_probs The two probabilities defining the fallback spread, used
#'   only when the interquartile range is zero.
#'
#' @return A numeric vector the same length as \code{y}, carrying \code{centre}
#'   and \code{scale} attributes for \code{\link{back_transform}}.
#'
#' @seealso \code{\link{back_transform}}, \code{\link{rre}}
#'
#' @examples
#' y <- c(stats::rlnorm(100), 5000)   # one extreme value
#' z <- robust_normalize(y)
#' c(centre = attr(z, "centre"), scale = attr(z, "scale"))
#'
#' @export
robust_normalize <- function(y, fallback_probs = c(0.10, 0.90)) {
  centre <- median(y, na.rm = TRUE)
  scale  <- IQR(y, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    q <- quantile(y, fallback_probs, na.rm = TRUE, names = FALSE)
    scale <- q[2] - q[1]
    warning("IQR was zero; fell back to the ",
            paste0(fallback_probs * 100, collapse = "-"), " percentile spread.")
  }
  if (!is.finite(scale) || scale <= 0)
    stop("Variable has no usable spread; check for a constant or all-censored column.")
  structure((y - centre) / scale, centre = centre, scale = scale)
}

#' @title Undo a robust rescaling
#'
#' @description
#' Returns values rescaled by \code{\link{robust_normalize}} to their original
#' units, for plotting fitted trajectories on the scale a reader recognises.
#'
#' @param y_scaled Numeric vector on the rescaled scale.
#' @param centre The \code{centre} attribute from \code{\link{robust_normalize}}.
#' @param scale The \code{scale} attribute from \code{\link{robust_normalize}}.
#'
#' @return A numeric vector in the original units.
#'
#' @seealso \code{\link{robust_normalize}}
#'
#' @examples
#' y <- stats::rlnorm(50)
#' z <- robust_normalize(y)
#' all.equal(y, back_transform(z, attr(z, "centre"), attr(z, "scale")),
#'           check.attributes = FALSE)
#'
#' @export
back_transform <- function(y_scaled, centre, scale) y_scaled * scale + centre

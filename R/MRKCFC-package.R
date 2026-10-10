#' @keywords internal
"_PACKAGE"

#' @importFrom stats ave median IQR quantile rchisq rnorm runif spline uniroot
NULL

# These are imported here rather than left bare because a bare `median` or
# `quantile` resolves at run time by searching the attached packages. That
# works only while `stats` is attached and nothing earlier on the search path
# defines the name; inside a parallel worker, where the environment is rebuilt,
# it may not resolve at all. robust_normalize() and robust_stage2() both depend
# on `median`, `IQR` and `quantile`, and both run inside the K-centers loop.
# Declaring the imports once here covers every call site, including new ones.

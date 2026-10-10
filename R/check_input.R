#' Validate multivariate curve inputs
#'
#' Checks that \code{Ly} and \code{Lt} have the shape \code{\link{MKCFC}}
#' expects and returns the number of subjects.
#'
#' @param Ly A list of \code{P} elements, one per variable, each a list of
#'   \code{n} numeric vectors.
#' @param Lt A list of the same shape holding observation times. Note this is a
#'   list of \code{P}, not of \code{n}.
#'
#' @return The number of subjects, invisibly to the caller's use. Stops with a
#'   message naming the violation if the shapes disagree.
#'
#' @seealso \code{\link{MKCFC}}
#' @export
check_input <- function(Ly, Lt){
  .mrkcfc_validate_multivariate(Ly, Lt)
}

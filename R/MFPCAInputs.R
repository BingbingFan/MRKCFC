#' Convert a long data frame to the list form MKCFC expects
#'
#' Reshapes one row per observation into the \code{Ly} and \code{Lt} lists used
#' throughout the package.
#'
#' @param id Name of the subject identifier column in \code{data}.
#' @param y Character vector naming the \code{P} outcome columns.
#' @param timevar Name of the observation-time column.
#' @param data A long-format data frame containing those columns.
#'
#' @return A list with \code{Ly}, \code{Lt} and \code{id}, the last giving the
#'   subject identifiers in the order the lists use.
#'
#' @seealso \code{\link{MKCFC}}, which calls this when given the data-frame
#'   input form.
#' @export
MFPCAInputs <- function(id, y, timevar, data){
  required <- c(id, y, timevar)
  if (!all(required %in% names(data))) {
    stop("id, y, and timevar must name columns in data.", call. = FALSE)
  }
  uid <- unique(data[[id]])
  Ly <- vector("list", length(y))
  Lt <- vector("list", length(y))
  for (j in seq_along(y)) {
    Ly[[j]] <- vector("list", length(uid))
    Lt[[j]] <- vector("list", length(uid))
    for (i in seq_along(uid)) {
      subject <- data[data[[id]] == uid[[i]], c(y[[j]], timevar), drop = FALSE]
      subject <- stats::na.omit(subject)
      subject <- subject[order(subject[[timevar]]), , drop = FALSE]
      if (nrow(subject) == 0L) {
        stop("Each subject must have at least one observed value for every outcome.", call. = FALSE)
      }
      Ly[[j]][[i]] <- subject[[y[[j]]]]
      Lt[[j]][[i]] <- subject[[timevar]]
    }
  }
  list(id = as.list(uid), Ly = Ly, Lt = Lt)
}

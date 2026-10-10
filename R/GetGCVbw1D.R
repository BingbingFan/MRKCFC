#' @title Bandwidth selection through GCV for one-dimension cases
#'
#' @description Bandwidth selection through generalized cross-validation (GCV) for one-dimension cases.
#'
#' @param Lt A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the observation time in ascending order for each subject.
#' @param Ly A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the measurements of each subject at the observation time correspond to \code{Lt}.
#' @param kern A \code{character} denoting the kernel type; 'epan'(Epanechnikov), 'unif'(Uniform), 'quar'(Quartic), 'gauss'(Gaussian).
#' @param dataType A \code{character} denoting the data type; 'Sparse'-default, 'Dense'.
#'
#' @return A scalar denoting the optimal bandwidth.
#' @export
#' @examples
#' # GenDataKL() is the data generator from the upstream KFPCA package. Only
#' # the KFPCA runtime is vendored here, so this example is shown rather than
#' # run; install KFPCA to execute it.
#' \dontrun{
#' # Generate data
#' n <- 100
#' interval <- c(0, 10)
#' lambda_1 <- 9 #the first eigenvalue
#' lambda_2 <- 1.5 #the second eigenvalue
#' eigfun <- list()
#' eigfun[[1]] <- function(x){cos(pi * x/10)/sqrt(5)}
#' eigfun[[2]] <- function(x){sin(pi * x/10)/sqrt(5)}
#' score <- cbind(rnorm(n, 0, sqrt(lambda_1)), rnorm(n, 0, sqrt(lambda_2)))
#' DataNew <- GenDataKL(n, interval = interval, sparse = 6:8, regular = FALSE,
#'                      meanfun = function(x){0}, score = score,
#'                      eigfun = eigfun, sd = sqrt(0.1))
#' # Optimal bandwidth for mean function estimate
#' bwOpt <- GetGCVbw1D(DataNew$Lt, DataNew$Ly, kern = "epan")
#' }
GetGCVbw1D <- function(Lt, Ly, kern, dataType = "Sparse"){

  return(GCVLwls1D1(Ly, Lt, kern, dataType = dataType)$bOpt)

}

#' @title Bandwidth selection through GCV for two-dimension cases
#'
#' @description Bandwidth selection through generalized cross-validation (GCV) for two-dimension cases.
#'
#' @param tPairs A \code{matrix} with two columns containing the pairs of time points.
#' @param yin A \code{vector} denoting the corresponding values.
#' @param Lt A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the observation time in ascending order for each subject.
#' @param kern A \code{character} denoting the kernel type; 'epan'(Epanechnikov), 'unif'(Uniform), 'quar'(Quartic), 'gauss'(Gaussian).
#' @param ObsGrid A \code{vector} containing all observation grids in ascending order.
#' @param RegGrid A \code{vector} of the equally spaced time points in the support interval.
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
#' # Optimal bandwidth for the estimate of
#' # E{X(s)X(t)} = cov(X(s), X(t)) + mu(s) * mu(t)
#' xin2D <- NULL
#' yin2D <- NULL
#' for(i in 1:n){
#'   xin2D <- rbind(xin2D, t(utils::combn(DataNew$Lt[[i]], 2)))
#'   yin2D <- rbind(yin2D, t(utils::combn(DataNew$Ly[[i]], 2)))
#' }
#' tPairs <- xin2D
#' yin <- yin2D[,1] * yin2D[, 2]
#' bwOpt <- GetGCVbw2D(tPairs = tPairs, yin = yin, Lt = DataNew$Lt,
#'                     kern = "epan", ObsGrid = sort(unique(unlist(DataNew$Lt))),
#'                     RegGrid = seq(interval[1], interval[2], length.out = 51))
#' }
GetGCVbw2D <- function(tPairs, yin, Lt, kern, ObsGrid, RegGrid, dataType = "Sparse"){

  rcov <- list()
  rcov$tPairs <- tPairs
  rcov$cxxn <- yin

  return(GCVLwls2DV2(ObsGrid, RegGrid, rcov = rcov, kern = kern, Lt = Lt, dataType = dataType)$h)
}

#' @title Least square estimates of functional principal component scores
#'
#' @description Least square estimates (LSE) of functional principal component scores.
#'
#' @param Lt A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the observation time in ascending order for each subject.
#' @param Ly A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the measurements of each subject at the observation time correspond to \code{Lt}.
#' @param kern A \code{character} denoting the kernel type; 'epan'(Epanechnikov), 'unif'(Uniform), 'quar'(Quartic), 'gauss'(Gaussian).
#' @param bw A scalar denoting the bandwidth for mean function estimate.
#' @param FPC_dis A \code{nRegGrid} by \code{nK} \code{matrix} containing the eigenfunction estimates at \code{RegGrid}, where \code{nRegGrid} is the length of \code{RegGrid} and \code{nK} is the number of FPCs.
#' @param RegGrid A \code{vector} of the equally spaced time points in the support interval.
#' @param more Logical; If \code{FALSE}, only the estimates of FPC scores are returned. If \code{TRUE}, the mean function estimates and the eigenfunction estimates at all observation time points are also returned.
#'
#' @return If \code{more = FALSE}, a \emph{n} by \code{nK} \code{matrix} containing the estimates of the FPC scores is returned, where \emph{n} is the sample size. If \code{more = TRUE}, a \code{list} containing the following components is returned:
#' \item{score}{a \emph{n} by \code{nK} \code{matrix} containing the estimates of the FPC scores.}
#' \item{meanest_fine}{Mean function estimates at all observation time points.}
#' \item{FPC_dis_fine}{Eigenfunction estimates at all observation time points.}
#' @export
#'
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
#' DataNew <- GenDataKL(n, interval = interval, sparse = 3:5, regular = FALSE,
#'                      meanfun = function(x){0}, score = score,
#'                      eigfun = eigfun, sd = sqrt(0.1))
#' basis <- fda::create.bspline.basis(interval, nbasis = 13, norder = 4,
#'                               breaks = seq(0, 10, length.out = 11))
#' Klist <- KFPCA(DataNew$Lt, DataNew$Ly, interval, nK = 2, bw = 1,
#'                nRegGrid = 51, fdParobj = basis)
#' # Just an example to explain the use of FPCscoreLSE().
#' # One can obtain FPC scores estimates for KFPCA method
#' # by KFPCA() directly. Note that FPCscoreLSE() can also be used
#' # to estimate FPC scores for methods except KFPCA.
#' scoreKFPCA <- FPCscoreLSE(DataNew$Lt, DataNew$Ly, kern = "epan",
#'                           bw = Klist$bwmean, FPC_dis = Klist$FPC_dis,
#'                           RegGrid = seq(interval[1], interval[2], length.out = 51))
#' }
FPCscoreLSE <- function(Lt, Ly, kern, bw, FPC_dis, RegGrid, more = FALSE){

  # n: the sample size
  n <- length(Lt)

  # meanest: mean function estimates at all observation grids
  # ObsGrid: all observation grids in ascending order
  meanestList <- MeanObsGrid(Lt, Ly, kern, bw)
  meanest <- meanestList$mean
  ObsGrid <- meanestList$ObsGrid

  # nK: the number of FPCs
  # nObsGrid: the number of observation grids
  # FPC_dis_aug: A nObsGrid by nK matrix containing eigenfunction estimates at various observation grids
  nK <- ncol(FPC_dis)
  nObsGrid <- length(ObsGrid)
  bwFPC <- NULL
  FPC_dis_aug <- matrix(0, ncol = nK, nrow = nObsGrid)
  for(i in 1:nK){
    bwFPC[i] <- GetGCVbw1D(Lt = list(RegGrid), Ly = list(FPC_dis[,i]), kern = "epan",
                           dataType = "Dense")
    FPC_dis_aug[,i] <- fdapace::Lwls1D(bwFPC[i], kernel_type = "epan", xin = RegGrid, yin = FPC_dis[,i],
                                       xout = ObsGrid)
  }

  # mean_ind: a list of vectors. The i-th vector contains the mean estimates at the observation
  # time of the i-th subject.
  # phi: a list of matrices. The i-th matrix contains the eigenfunction estimates at the
  # observation time of the i-th subject.
  # score: a n by nK matrix containing the estimates of the FPC scores. The (i, j)-th element is
  # the j-th FPC score estimate of the i-th subject.
  mean_ind <- list()
  phi <- list()
  score <- matrix(0, ncol = nK, nrow = n)
  for(i in 1:n){
    id <- sapply(Lt[[i]], function(x){which(ObsGrid == x)})
    mean_ind[[i]] <- meanest[id]
    phi[[i]] <- FPC_dis_aug[id,]
    response <- Ly[[i]] - mean_ind[[i]]
    score[i,] <- as.vector(.mrkcfc_solve_covariance(
      crossprod(phi[[i]]),
      crossprod(phi[[i]], response)
    ))
  }

  if(more == FALSE){
    return(score)
  }else{
    ret <- list()
    ret$score <- score
    ret$meanest_fine <- meanest
    ret$FPC_dis_fine <- FPC_dis_aug
    return(ret)
  }

}



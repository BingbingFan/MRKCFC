#' @title Kendall Functional Principal Component Analysis (KFPCA) for sparse design
#'
#' @description KFPCA for non-Gaussian functional data with sparse design or longitudinal data.
#'
#' @param Lt A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the observation time in ascending order for each subject.
#' @param Ly A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the measurements of each subject at the observation time correspond to \code{Lt}.
#' @param interval A \code{vector} of length two denoting the supporting interval.
#' @param dataType A \code{character} denoting the data type; 'Sparse'-default, 'Dense'.
#' @param nK An integer denoting the number of FPCs.
#' @param kern A \code{character} denoting the kernel type for the Nadaraya-Watson estimators; 'epan'(Epanechnikov)-default, 'unif'(Uniform), 'quar'(Quartic), 'gauss'(Gaussian).
#' @param bw A scalar denoting the bandwidth for the Nadaraya-Watson estimators.
#' @param kernK A \code{character} denoting the kernel type for the estimation of the Kendall's tau function; 'epan'(Epanechnikov)-default, 'unif'(Uniform), 'quar'(Quartic), 'gauss'(Gaussian).
#' @param bwK The bandwidth for the estimation of the Kendall's tau function. If \code{is.numeric(bwK) == T}, \code{bwK} is exactly the bandwidth. If \code{bwK == "GCV"}, the bandwidth is chosen by GCV. (default: "GCV")
#' @param kernmean A \code{character} denoting the kernel type for the estimation of the mean function; 'epan'(Epanechnikov)-default, 'unif'(Uniform), 'quar'(Quartic), 'gauss'(Gaussian).
#' @param bwmean The bandwidth for the estimation of the mean function. If \code{is.numeric(bwmean) == T}, \code{bwmean} is exactly the bandwidth. If \code{bwmean == "GCV"}, the bandwidth is chosen by GCV. (default: "GCV")
#' @param nRegGrid An integer denoting the number of equally spaced time points in the supporting interval. The eigenfunctions and mean function are estimated at these equally spaced time points.
#' @param fdParobj A functional parameter object for the smoothing of the eigenfunctions. For more detail, see \code{\link[fda]{smooth.basis}}.
#' @param more Logical; If \code{FALSE}, estimates of FPC scores and predictions of trajectories are not returned.
#'
#' @return A \code{list} containing the following components:
#' \item{ObsGrid}{A \code{vector} containing all observation time points in ascending order.}
#' \item{RegGrid}{A \code{vector} of the equally spaced time points in the support interval.}
#' \item{bwmean}{A scalar denoting the bandwidth for the mean function estimate.}
#' \item{kernmean}{A \code{character} denoting the kernel type for the estimation of the mean function}
#' \item{bwK}{A scalar denoting the bandwidth for the Kendall's tau function estimate.}
#' \item{kernK}{A \code{character} denoting the kernel type for the estimation of the Kendall's tau function}
#' \item{mean}{A \code{vector} of length \code{nRegGrid} denoting the mean function estimate.}
#' \item{KendFun}{A \code{nRegGrid} by \code{nRegGrid} \code{matrix} denoting the Kendall's tau function estimate.}
#' \item{FPC_dis}{A \code{nRegGrid} by \code{nK} \code{matrix} containing the eigenfunction estimates at \code{RegGrid}.}
#' \item{FPC_smooth}{A functional data object for the eigenfunction estimates.}
#' \item{score}{A \emph{n} by \code{nK} \code{matrix} containing the estimates of the FPC scores, where \emph{n} is the sample size. The results are returned when \code{more = TRUE}.}
#' \item{X_fd}{A functional data object for the prediction of trajectories. The results are returned when \code{more = TRUE}.}
#' \item{Xest_ind}{A \code{list} containing the prediction of each trajectory at their own observation time points. The results are returned when \code{more = TRUE}.}
#' \item{Lt}{The input 'Lt'.}
#' \item{Ly}{The input 'Ly'.}
#' \item{CompTime}{A scalar denoting the computation time.}
#' @export
#'
#' @references
#' \cite{Rou Zhong, Shishi Liu, Haocheng Li, Jingxiao Zhang (2021). "Robust Functional Principal Component Analysis for Non-Gaussian Longitudinal Data." Journal of Multivariate Analysis, https://doi.org/10.1016/j.jmva.2021.104864.}
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
#' DataNew <- GenDataKL(n, interval = interval, sparse = 6:8, regular = FALSE,
#'                      meanfun = function(x){0}, score = score,
#'                      eigfun = eigfun, sd = sqrt(0.1))
#' basis <- fda::create.bspline.basis(interval, nbasis = 13, norder = 4,
#'                               breaks = seq(0, 10, length.out = 11))
#' # KFPCA
#' Klist <- KFPCA(DataNew$Lt, DataNew$Ly, interval, nK = 2, bw = 1,
#'                nRegGrid = 51, fdParobj = basis)
#' plot(Klist$FPC_smooth)
#' }
KFPCA <- function(Lt, Ly, interval, dataType = "Sparse", nK,
                  kern = "epan", bw, kernK = "epan", bwK = "GCV",
                  kernmean = "epan", bwmean = "GCV",
                  nRegGrid, fdParobj, more = TRUE){
  startFPCA <- Sys.time() #Starting time for KFPCA

  # n: sample size
  n <- length(Lt)

  # Xlist: the list of n matrices, where n is the sample size. The i-th matrix contains
  # the observations of the i-th subject and the Nadaraya-Watson estimators of the other
  # subjects at the observation time for the i-th subject
  Xlist <- GetMatList(Lt, Ly, kern, bw)

  # For each matrix in Xlist, drop out the rows that contain NaNs
  naid <- which(unlist(lapply(Xlist, function(x){length(which(is.nan(x) == T))})) != 0)
  for(i in naid){
    narow <- which(apply(Xlist[[i]], 1, function(x){length(which(is.nan(x) == T))}) != 0)
    Xlist[[i]] <- Xlist[[i]][-narow,]
  }
  sub_id <- which(unlist(lapply(Xlist, is.matrix)) == T) 
  #the subjects that can be used to compute the raw Kendall's tau covariances

  # Lk: the list of n matrices, where n is the sample size. The i-th matrix is the raw Kendall's
  # tau covariance obtained from the i-th subject
  Lk <- NULL
  for(i in sub_id){
    mi <- ncol(Xlist[[i]])
    ni <- nrow(Xlist[[i]])
    Lk[[i]] <- matrix(0, ncol = mi, nrow = mi)
    for(j in (2:ni)){
      XX <- Xlist[[i]][1,] - Xlist[[i]][j,]
      Lk[[i]] <- Lk[[i]] + (XX %*% t(XX))/(normq(XX, 2))^2 * length(XX)
    }
    Lk[[i]] <- Lk[[i]]/(ni - 1)
  }

  # xintotal: a matrix of the assembled time pairs for the raw Kendall's tau covariance
  xin <- NULL
  for(i in sub_id){
    xin <- rbind(xin, cbind(i, t(utils::combn(Lt[[i]], 2))))
  }
  xinsym <- cbind(xin[,1], xin[,3], xin[,2])
  xintotal <- rbind(xin, xinsym)

  # yin: the raw Kendall's tau covariance for various assembled time pairs
  yin <- NULL
  for(i in sub_id){
    y_id <- which(lower.tri(Lk[[i]]) == T)
    yin <- c(yin, Lk[[i]][y_id])
  }
  yin <- c(yin, yin)

  # Bandwidth for Kendall's tau function estimate
  gridobs <- sort(unique(c(unlist(Lt))))
  gridreg <- seq(interval[1], interval[2], length.out = nRegGrid)
  if(is.numeric(bwK)){
    bw2 <- bwK
  }else if(bwK == "GCV"){
    bw2 <- GetGCVbw2D(xintotal[,2:3], yin, Lt, kernK, gridobs, gridreg, dataType = dataType)
  }else{
    stop('bwK must be a scalar or the character denoting GCV !')
  }

  # Local linear estimate of the Kendall's tau function
  Kend <- fdapace::Lwls2D(bw = bw2, xin = xintotal[,2:3], yin = yin,
                          xout1 = gridreg, xout2 = gridreg)

  # FPC_dis: the first nK eigenfunction estimates at regular grids
  # FPC_smooth: the functional data object for the first nK eigenfunction estimates
  #
  # The eigenvalues were previously computed here and DISCARDED, which left no
  # way to compute a fraction of variance explained on the robust branch and so
  # no way to choose nK from the data -- fdapace returns $lambda and $cumFVE,
  # KFPCA returned neither. They are kept below as $values.
  #
  # NOTE for anyone using them: these are eigenvalues of the DISCRETISED
  # covariance surface, so they carry a grid factor. That factor cancels in any
  # ratio, which is all FVE needs. A smoothed PACE-type surface is also not
  # guaranteed positive semi-definite, so NEGATIVE eigenvalues occur and must be
  # truncated at zero before any cumulative share is taken.
  .eig <- eigen(Kend)
  FPC_dis <- apply(as.matrix(.eig$vectors[,1:nK]), 2, unitize, grid = diff(gridreg)[1])
  FPC_smooth <- fda::smooth.basis(gridreg, FPC_dis, fdParobj)$fd

  # Bandwidth for mean function estimates
  if(is.numeric(bwmean)){
    bw3 <- bwmean
  }else if(bwmean == "GCV"){
    bw3 <- GetGCVbw1D(Lt, Ly, kernmean, dataType = dataType)
  }else{
    stop('bwmean must be a scalar or the character denoting GCV !')
  }
  meanest <- MeanEst(Lt, Ly, kern = kernmean, bw3, gridreg)$mean

  # The return object
  ret <- list()
  ret$ObsGrid <- gridobs
  ret$RegGrid <- gridreg
  ret$bwmean <- bw3
  ret$kernmean <- kernmean
  ret$bwK <- bw2
  ret$kernK <- kernK
  ret$mean <- meanest
  ret$KendFun <- Kend
  ret$FPC_dis <- FPC_dis
  ret$values <- .eig$values      # full spectrum; see the note above
  ret$FPC_smooth <- FPC_smooth
  ret$Lt <- Lt
  ret$Ly <- Ly

  class(ret) <- 'KFPCA'

  if(more){
    # FPC scores
    score_ret <- FPCscoreLSE(Lt, Ly, kern = kernmean, bw3, FPC_dis, gridreg, more = TRUE)
    FPC_score <- score_ret$score

    # X_dis: the prediction of trajectories at regular grids
    # X_fd: the functional data object for the prediction of trajectories
    X_dis <- FPC_dis %*% t(FPC_score) + meanest
    X_fd <- fda::smooth.basis(gridreg, X_dis, fdParobj)$fd

    # X_dis_fine: the prediction of trajectories at all observation time points
    # Xest_ind: A list containing the prediction of each trajectory at their own
    # observation time points
    X_dis_fine <- score_ret$FPC_dis_fine %*% t(FPC_score) + score_ret$meanest_fine
    Xest_ind <- list()
    for(i in 1:n){
      Xest_ind[[i]] <- X_dis_fine[match(Lt[[i]], gridobs), i]
    }

    #Computation time
    endFPCA <- Sys.time()
    CompTime <- endFPCA - startFPCA

    ret$score <- FPC_score
    ret$X_fd <- X_fd
    ret$Xest_ind <- Xest_ind
    ret$CompTime <- CompTime

    return(ret)
  }else{
    #Computation time
    endFPCA <- Sys.time()
    CompTime <- endFPCA - startFPCA

    ret$CompTime <- CompTime

    return(ret)
  }

}



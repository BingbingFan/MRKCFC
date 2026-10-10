#' Modified multivariate silhouette for a functional partition
#'
#' A silhouette adapted to sparsely observed curves. For each subject, the
#' distance to a class is the mean squared \eqn{L_2} deviation of its observed
#' values from that class's estimated mean function, evaluated at the subject's
#' own observation times. \eqn{a_i} is that distance to the subject's own class,
#' \eqn{b_i} the smallest to any other, and the score averages
#' \eqn{(b_i - a_i) / \max(a_i, b_i)}.
#'
#' @section What it can and cannot answer: the criterion is defined only for
#' \strong{two or more classes} -- \eqn{b_i} requires a nearest other class --
#' so it cannot report the \emph{absence} of clustering. It answers how many
#' classes there are given that there are at least two. A low value is
#' consistent with weakly separated genuine classes and equally consistent with
#' no class structure at all, and this criterion cannot distinguish those. Where
#' that distinction matters, rest the conclusion on an external criterion rather
#' than on internal fit.
#'
#' Because distances are taken to class \emph{mean functions}, the criterion
#' rewards separation in level and penalises separation in shape or in rate of
#' change. Two partitions that differ in which of those they capture are not
#' ranked neutrally by it.
#'
#' @param data A long-format data frame with one row per observation, holding
#'   the identifier, value, time and class columns named below, the class column
#'   being \code{g}.
#' @param class Name of the class column in \code{data}.
#' @param id Name of the subject identifier column.
#' @param y Name of the observed value column.
#' @param time Name of the observation time column.
#' @param robust Logical, whether class mean functions are estimated robustly.
#' @param option Optional list of estimation options passed to the mean-function
#'   estimator.
#'
#' @return A list with components
#'   \item{a_i}{per-subject mean squared deviation from its own class mean.}
#'   \item{b_i}{per-subject mean squared deviation from the nearest other class
#'     mean.}
#'   \item{silhouette}{the average of \eqn{(b_i - a_i) / \max(a_i, b_i)} over
#'     subjects. Higher is better.}
#'
#' @seealso \code{\link{MKCFC}} for the fit whose partition this scores.
#'
#' @export
silhouette_score <- function(data, class, id, y, time,  robust = FALSE, option = NULL ){
  # require(fdapace)
  # if(robust){
  #   require(KFPCA)
  # }
  # require(dplyr)
  data <- data[,c( id,class, time, y)]
  names(data) <- c( "id","g", "time", "y")
  K <- length(unique(data$g))
  n <- length(unique(data$id))
  # FPCA data
  dtk <- lapply(1:K, function(x){
    # Explicit indexing, not subset(). `subset()` resolves `g` by non-standard
    # evaluation, which R CMD check cannot see -- and if the column were ever
    # renamed it would look outside the frame rather than erroring.
    data[data$g == x, , drop = FALSE]
  })

  fdata_k <- lapply(1:K, function(x) fdapace::MakeFPCAInputs(IDs = dtk[[x]]$id, yVec = dtk[[x]]$y, tVec = dtk[[x]]$time ))
  # estimate mean function using FPCA or KFPCA
  mu_func <- function(dt, option){
    if(robust){
      # bwK was previously omitted here, and ONLY here: MFPCA_k(),
      # clust_update() and CH_sparse() all pass optns$bwK. Under robust = TRUE
      # that silently smoothed the covariance with KFPCA's default "GCV"
      # instead of the configured bandwidth, so the silhouette's mean functions
      # were estimated under different smoothing from the model being scored.
      FPCA_obj <- KFPCA(Lt = dt$Lt,  Ly = dt$Ly, interval = option$interval, nK = option$nK, bw = option$bw,
                                bwK = option$bwK,
                                nRegGrid = option$nRegGrid, fdParobj = option$basis)
      mu <- FPCA_obj$mean
      Grid <- FPCA_obj$RegGrid
    }else{
      FPCA_obj <- fdapace::FPCA(Ly = dt$Ly, Lt = dt$Lt, optns = option)
      mu <- FPCA_obj$mu
      Grid <- FPCA_obj$workGrid
    }
    return(list(mu = mu, Grid = Grid))
  }

  FPCA_res <- lapply(1:K, function(x) mu_func(dt = fdata_k[[x]], option = option) )
  mus <- lapply(1:K, function(x) FPCA_res[[x]]$mu )
  Grids <- lapply(1:K, function(x) FPCA_res[[x]]$Grid )
  # within-cluster mean sum of square
  # for subjects i , calculate L2 distance of i with u_k(t)
  fdata <- fdapace::MakeFPCAInputs(IDs = data$id, yVec = data$y, tVec = data$time )

  L2_dis <- function(Lyi, Lti, mus, Grids){
    mean_i_k <-  lapply(1:K, function(x) spline(Grids[[x]] ,mus[[x]], xout = Lti)$y)
    d_L2 <- lapply(1:K, function(x) (Lyi - mean_i_k[[x]])^2 )
    d_L2 <- unlist(lapply(1:K, function(x) mean(d_L2[[x]])))
  }

  L2_i <- lapply(1:n, function(x) L2_dis(Lyi = fdata$Ly[[x]], Lti = fdata$Lt[[x]], mus, Grids) )
  #   cal calculate b_i and a_i
  uid <- data.frame(id = unlist(fdata$Lid))
  class_i <- data.frame(unique(data[, c("id","g")]))
  class_i <- merge(uid, class_i,  by="id", all.x = T )
  # bi ai fun
  a_i <- unlist(lapply(1:n, function(x) L2_i[[x]][ class_i$g[x] ]  ))
  b_i <- unlist(lapply(1:n, function(x) min(L2_i[[x]][ - class_i$g[x] ])  ))

  silhouette <- mean( unlist( lapply(1:n, function(x) (b_i[x] - a_i[x])/max(b_i[x], a_i[x]) ) ) )

  return(list(a_i = a_i, b_i = b_i , silhouette = silhouette))
}

#' Reconstruct one subject's trajectories under a fitted class
#'
#' Given a class's multivariate FPCA fit, reconstructs a single subject's
#' curves at that subject's own observation times. This is the quantity the
#' K-centres loop compares across classes when reassigning: a subject is
#' assigned to the class whose reconstruction fits it best.
#'
#' The step is a \emph{reconstruction}, not a forecast. It fits the observed
#' trajectory under a class-specific model; it does not predict beyond the
#' observation window.
#'
#' @section Known limitation: fails on a subject with only one observation on
#' any variable, because the per-variable score matrix collapses to a bare
#' vector rather than a one-row matrix. Apply the model where every subject has
#' at least two observations per variable; \code{\link{MKCFC}} is used that way
#' throughout.
#'
#' @param FPCA_k A list of \code{P} marginal FPCA fits for one class, as carried
#'   in \code{ref_hatyc} by \code{\link{MFPCA_k}}.
#' @param ith_yy A list of \code{P} numeric vectors: the subject's observed
#'   values on each variable.
#' @param target_t A list of \code{P} numeric vectors of the corresponding
#'   observation times.
#' @param robust Logical, whether \code{FPCA_k} came from the robust branch.
#'   Must match how the fits were produced, since the two carry differently
#'   named components.
#' @param optns Options for the marginal step; see \code{\link{MKCFC}}.
#' @param mfpca_scatter Scatter for the second stage, \code{"sign"} or
#'   \code{"moment"}.
#'
#' @return A list with components
#'   \item{fittedX}{a list of \code{P} numeric vectors, the reconstructed values
#'     at \code{target_t}.}
#'   \item{multi_score}{the subject's multivariate scores under this class.}
#'   \item{multi_eigen}{the multivariate eigenfunctions used.}
#'
#' @seealso \code{\link{MFPCA_k}} for the fits this consumes, and
#'   \code{\link{MKCFC}} for the loop that calls it.
#'
#' @export
predicted_yhat <- function(FPCA_k, ith_yy,  target_t, robust, optns,
                           mfpca_scatter = c("sign", "moment") ){
  # suppressMessages(require(KFPCA))
  J <- length(FPCA_k)
  ## i th yy : curve that I want to predict
  if(robust){

    n_M <- lapply(1:length(ith_yy), function(x) dim(FPCA_k[[x]]$score)[2] )
    phi <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$FPC_smooth )
    mu <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$mean )
    workGrid <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$RegGrid )

    prd <- lapply(seq_len(J), function(x) stats::predict(
      FPCA_k[[x]],
      newLt = list(target_t[[x]]),
      newLy = list(ith_yy[[x]]),
      nK = optns$nK,
      more = TRUE
    ))
    xi_Est_i <- lapply(1:J, function(x) prd[[x]]$score_new)
    pi_temp <- lapply(1:J, function(x) prd[[x]]$FPC_dis_new)
    # mfpca
    xiEst <- lapply(1:J, function(x) FPCA_k[[x]]$score)

    full_xiEst <- do.call(cbind, xiEst)
    rotation <- .mrkcfc_mfpca_rotation(full_xiEst, mfpca_scatter)
    mu_obs <- lapply(seq_len(J), function(j)
      spline(workGrid[[j]], mu[[j]], xout = target_t[[j]])$y)
    recon <- .mrkcfc_mfpca_reconstruct(rotation, pi_temp, xi_Est_i, n_M, mu_obs)
    fittedX <- recon$fittedX
    multi_score <- recon$multi_score
    multi_eigen <- recon$multi_eigen

  }else{
    n_M <- lapply(1:length(ith_yy), function(x) length(FPCA_k[[x]]$lambda) )
    lambda <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$lambda )
    phi <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$phi )
    mu <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$mu )
    workGrid <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$workGrid )
    sig2 <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$sigma2 )

    xiEst <- fittedX <- xi_Est_i <-  list()

    pi_observed <- vector("list", J)
    for(j in seq_len(J)){
      pi_observed[[j]] <- vapply(seq_len(n_M[[j]]), function(l) {
        stats::spline(workGrid[[j]], phi[[j]][, l], xout = target_t[[j]])$y
      }, numeric(length(ith_yy[[j]])))
      mu_ij <- stats::spline(workGrid[[j]], mu[[j]], xout = target_t[[j]])$y
      xi_Est_i[[j]] <- .mrkcfc_nonrobust_scores(
        pi_observed[[j]],
        lambda[[j]],
        sig2[[j]],
        ith_yy[[j]] - mu_ij
      )
      xiEst[[j]]  <-  FPCA_k[[j]]$xiEst
    }

    full_xiEst <- do.call(cbind, xiEst)
    rotation <- .mrkcfc_mfpca_rotation(full_xiEst, mfpca_scatter)
    mu_obs <- lapply(seq_len(J), function(j)
      stats::spline(workGrid[[j]], mu[[j]], xout = target_t[[j]])$y)
    recon <- .mrkcfc_mfpca_reconstruct(rotation, pi_observed, xi_Est_i, n_M, mu_obs)
    fittedX <- recon$fittedX
    multi_score <- recon$multi_score
    multi_eigen <- recon$multi_eigen

  }

  return(list(fittedX=fittedX, multi_score = multi_score, multi_eigen = multi_eigen ))
}

#' Multivariate functional principal component analysis within one class
#'
#' Performs the two-stage multivariate FPCA of Happ and Greven (2018) on the
#' subjects assigned to a single class: a marginal functional principal
#' component analysis of each variable, then an eigendecomposition of the
#' stacked univariate scores.
#'
#' @section The two stages are controlled separately: \code{robust} selects the
#' marginal estimator, Kendall FPCA when \code{TRUE} and ordinary FPCA when
#' \code{FALSE}. \code{mfpca_scatter} selects the scatter whose eigenvectors
#' give the second-stage rotation. The two are independent, and the mixed
#' setting \code{robust = TRUE} with \code{mfpca_scatter = "moment"} recovers
#' the eigenfunction subspace worse than either consistent choice: the robust
#' marginal step leaves non-negligible variance in components the non-robust one
#' drives to nearly zero, and a moment scatter weights those while a sign
#' scatter does not.
#'
#' @section On using \code{multi_rotation}: the multivariate eigenfunctions are
#' obtained by stacking the elements, \eqn{\psi_m^{(j)} = \Phi_j C_{[j], m}} in
#' the notation of Happ and Greven's Equation (9), where \eqn{C} is
#' \code{multi_rotation} and \eqn{[j]} is variable \code{j}'s block of rows.
#' \strong{The elements must be stacked, not scored one at a time.} Taken
#' variable by variable, \eqn{\Phi_j C_{[j], m}} spans the same subspace as
#' \eqn{\Phi_j} itself, so any comparison made that way is blind to the rotation
#' -- which is the entire content of the second stage.
#'
#' Block sizes are \emph{not} guaranteed equal: with \code{robust = FALSE},
#' \pkg{fdapace} may select fewer components for one variable than another. Read
#' the per-variable component count from the fits in \code{ref_hatyc} rather
#' than assuming \code{optns$nK} for every block.
#'
#' @param u Index of the class to fit, in \code{1:K}.
#' @param Ly A list of \code{P} elements, one per variable, each a list of
#'   \code{n} numeric vectors of observed values.
#' @param Lt A list of the same shape holding observation times. Note this is a
#'   list of \code{P}, not of \code{n}.
#' @param clust_f Integer vector of length \code{n} giving each subject's class.
#'   Only subjects with \code{clust_f == u} are used.
#' @param optns Options for the marginal step. The two branches take different
#'   lists: see \code{\link{MKCFC}}.
#' @param robust Logical, whether the marginal step uses Kendall FPCA.
#' @param mfpca_scatter Scatter for the second stage, \code{"sign"} for the
#'   spatial-sign covariance matrix or \code{"moment"} for the sample second
#'   moment.
#'
#' @return A list with components
#'   \item{multi_score}{the multivariate scores of this class's subjects, the
#'     stacked univariate scores post-multiplied by \code{multi_rotation}.}
#'   \item{multi_rotation}{the rotation \eqn{C}, of dimension
#'     \code{sum(block sizes)} square, whose columns give the multivariate
#'     eigenfunctions when applied blockwise to the marginal eigenfunctions.}
#'   \item{ref_hatyc}{the \code{P} marginal fits, from \code{KFPCA} when
#'     \code{robust = TRUE} and \code{fdapace::FPCA} otherwise. Carries the
#'     eigenfunctions, means and per-variable component counts.}
#'
#' @references
#' Happ, C. and Greven, S. (2018). Multivariate functional principal component
#' analysis for data observed on different (dimensional) domains. \emph{Journal
#' of the American Statistical Association}, 113, 649--659.
#'
#' @seealso \code{\link{MKCFC}}, which calls this once per class each iteration.
#'
#' @export
MFPCA_k <- function(u, Ly, Lt,  clust_f, optns, robust,
                    mfpca_scatter = c("sign", "moment") ){
  xiEst <- multi_score <- NULL

  ClustIds <- .mrkcfc_cluster_ids(clust_f, length(unique(clust_f)))
  J <- length(Ly)
  ref_hatyc <- list()
  uu <- ClustIds[[u]]
  yy <- list()
  tt <- list()
  for(j in 1:J){
    yy[[j]] <-  lapply( uu,  function(x) Ly[[j]][[x]])
    tt[[j]] <- lapply( uu,  function(x) Lt[[j]][[x]])
  }

  if(robust){
    ref_hatyc <- lapply(1:J, function(x)
      KFPCA(Lt = tt[[x]],  Ly = yy[[x]], interval = optns$interval,
                    nK = optns$nK, bw = optns$bw, bwK = optns$bwK,
                    nRegGrid = optns$nRegGrid, fdParobj = optns$basis) )

    xiEst <- lapply(1:J, function(x) ref_hatyc[[x]]$score )
    n_M <- lapply(1:J, function(x) dim(ref_hatyc[[x]]$score)[2] )
    phi <- lapply(1:J, function(x) ref_hatyc[[x]]$FPC_dis )
    mux <- lapply(1:J, function(x) ref_hatyc[[x]]$mean )

  }else{

    ref_hatyc <- lapply(1:J, function(x)
      fdapace::FPCA(yy[[x]], tt[[x]], optns =  optns ))

    xiEst <- lapply(1:J, function(x) ref_hatyc[[x]]$xiEst )
    n_M <- lapply(1:J, function(x) ref_hatyc[[x]]$selectK )
    phi <- lapply(1:J, function(x) ref_hatyc[[x]]$phi )
    mux <- lapply(1:J, function(x) ref_hatyc[[x]]$mu )

  }
  full_xiEst <- do.call(cbind, xiEst)
  multi_rotation <- .mrkcfc_mfpca_rotation(full_xiEst, mfpca_scatter)
  multi_score <- full_xiEst %*% multi_rotation

  return(list(multi_score = multi_score,
              multi_rotation = multi_rotation,
              ref_hatyc = ref_hatyc  ))
}

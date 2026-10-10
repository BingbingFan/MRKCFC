# Bandwidth selection through generalized cross-validation (GCV) for two-dimension cases.
# The codes in this Rscript refer to the source code of fdapace package with a few modifications.

GCVLwls2DV2 <- function(obsGrid, regGrid, ngrid=NULL, dataType='Sparse', error=rcov$error, kern,
                        rcov, h0=NULL, verbose=FALSE, CV=FALSE, useBW1SE = FALSE, Lt) {

  # Returns: a list of length 2, containing the optimal bandwidth and the gcv score.
  # obsGrid: observation points. sort(unique(c(unlist(Lt))))

  ## Get minimal bandwidth and range
  r <- diff(range(obsGrid)) * sqrt(2) # sqrt(2) because the window is circular.
  minBW <- GetMinb(Lt, obsGrid = obsGrid, dataType = dataType)
  h0 <- minBW

  if (is.null(h0)){
    stop('the data is too sparse, no suitable bandwidth can be found! Try Gaussian Kernel instead!\n')
  }

  if (kern == 'gauss') {
    h0 = h0 * 0.2;
  }

  ## Get Candidate Bandwidths
  h0 <- min(h0, r/4)
  if (h0 < r/4) {
    q <- (r / (4 * h0)) ^ (1/9)
  } else if (h0 < r/2) {
    q <- (r / (2 * h0)) ^ (1/9)
  } else if (h0 < r) {
    q <- (r / h0) ^ (1/9)
  } else {
    stop('Data is too sparse. The minimal bandwidth is the range of data')
  }
  bw <- (q ^ seq(0,9,length.out = 10)) * h0 # from h0 to r / 4


  ## Set GCV/CV options
  opth <- h0

  leave <- FALSE
  iter <- 0
  maxIter <- 1

  minBWInvalid <- FALSE
  while (!leave && iter < maxIter) {
    if (minBWInvalid){
      minBW <- bw[1]
    }

    Scores <- matrix(Inf, nrow = length(bw), ncol = 2); colnames(Scores) <- c('SUM','SD');
    # try the bandwidths large to small in order to save time due to sparseness in the windows.
    for (i in rev(seq_along(bw))) {
      h <- bw[i]

      Scores[i,'SUM'] <- getGCVscoresV2(h, kern, rcov$tPairs, rcov$cxxn, regGrid=regGrid, verbose=verbose)

      if (is.infinite(Scores[i,'SUM'])) {
        minBWInvalid <- TRUE
        if (i < length(bw)) {
          if (minBWInvalid) {
            minBW <- bw[i + 1]
            minBWInvalid <- FALSE
          }
        }
        break; # This will help break out of the loop if the BW is too small to make sense
      }
    }

    if(is.infinite(min(Scores[,'SUM']))){
      opth <- max(bw)
      optgcv <- Inf
    } else {
      ind <- which.min(Scores[,'SUM'])
      opth <- bw[ind]
      optgcv <- Scores[ind,'SUM']
    }


    ## Check that what we found is coherent.
    if (opth >= r - 1e-12) {
      minBW <- r
      leave <- TRUE
      stop('Data is too sparse. The optimal bandwidth equals to the range of input time points. Try Gaussian kernel.')
    }
    if ( (abs(opth - max(bw)) > 1e-12) && !is.infinite(optgcv))
      leave <- TRUE
    else if (is.infinite(optgcv)) {
      if (verbose)
        warning('Data is too sparse, retry with larger bandwidths!')
      h0 <- max(bw) * 1.01
    } else if ( (abs(opth - max(bw)) > 1e-12) ) {
      warning('Optimal bandwidth not found in the candidate bandwidths. Retry with larger bandwidths')
      h0 <- max(bw)
    }
    if (!leave) {
      iter <- iter + 1
    }
  } # The "while (!leave && iter < maxIter) ..." end here

  ret <- list(h=opth, gcv=optgcv, minBW=minBW)
  if (CV != FALSE)
    names(ret)[2] <- 'cv'

  return(ret)

}

GetMinb <- function(t, obsGrid, dataType='Sparse', npoly=1, minUniqPts=3, minPts=6) {

  if (dataType == 'Sparse') {
    dstar <- Minb(obsGrid, 2 + npoly) # rough 1D initial value
    n_obs = length(obsGrid);
    tmp1 =  matrix( rep(0, n_obs^2), ncol = n_obs)

    # Find the pair against which we have measurements in the same curve
    for (i in 1:length(t)){
      idx = match( t[[i]], obsGrid)
      tmp1[idx, idx] = 1
    }
    res = tmp1 - diag(n_obs);
    ids = matrix(res > 0);
    b = matrix( rep(obsGrid, n_obs), nrow=n_obs);
    # Use half of the largest difference between two consequative points in the same
    # as curve as your candidate bandwith. We do no worry about the difference
    # between to [t_j(end) - t_{1+j}(1)] because this will be negative. This bandwidth tends to be conservative (too large).
    # dstar = max(dstar, max(diff(b[ids])/2)); # Original code
    dstar = max(dstar, stats::quantile( diff(b[ids]), 0.95)/2 ); # Fix to avoid outliers
  }else{
    dstar = Minb(obsGrid, 2 + npoly) * 1.5;
  }
  return(dstar)

}


getGCVscoresV2 <- function(bw, kern, xin, yin, win=NULL, regGrid, RSS=NULL, verbose=FALSE) {
  # RSS: for implementing GCV of binned rcov.
  # browser()

  if (is.null(win))
    win <- rep(1, length(yin))

  fit <- tryCatch(fdapace::Lwls2D(bw, kern, xin=xin, yin=yin, win=win, xout1=regGrid, xout2=regGrid), error=function(err) {
    if (verbose) {
      message('Invalid bandwidth. Try enlarging the window size.\n')
    }
    return(Inf)
  })

  # Catch
  if (is.infinite(fit[1]))
    return(Inf)

  # workaround for degenerate case.
  if (any(is.nan(fit)))
    return(Inf)

  obsFit <- pracma::interp2(regGrid, regGrid, fit, xin[, 1], xin[, 2])

  # residual sum of squares
  res <- sum((yin - obsFit) ^ 2 * win)
  if (!is.null(RSS))
    res <- res + sum(RSS)

  # kernel at x=0
  k0 <- KernelAt0(kern)
  N <- sum(win)
  r <- diff(range(xin[, 1]))
  bottom <- max(1 - 3 * (1 / N) * (r * k0 / bw)^2, 0)
  GCV <- res / bottom^2

  return(GCV)

}


KernelAt0 <- function(kern) {
  if (kern == 'quar')
    k0 <- 0.9375
  else if (kern == 'epan')
    k0 <- 0.75
  else if (kern == 'rect')
    k0 <- 0.5
  else if (kern == 'gausvar')
    k0 <- 0.498677850501791
  else if (kern == 'gauss')
    k0 <- 0.398942280401433
  else
    stop('Unknown kernel')

  return(k0)
}



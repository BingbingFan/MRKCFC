# Bandwidth selection through generalized cross-validation (GCV) for one-dimension cases.
# The codes in this Rscript refer to the source code of fdapace package with a few modifications.

GCVLwls1D1 <- function(Ly,Lt, kernel, npoly = 1, nder = 0, dataType = 'Sparse') {

  t = unlist(Lt);
  y = unlist(Ly)[order(t)];

  t = sort(t);

  # r = diff(range(t))
  N = length(t);
  r = t[N] - t[1];

  # Specify the starting bandwidth candidates
  if ( dataType == "Sparse") {
    dstar = Minb(t, npoly+2);
    if ( dstar > r*0.25){
      dstar = dstar * 0.75;
      warning( c( "The min bandwidth choice is too big, reduce to ", dstar, "!\n"))
    }
    h0 = 2.5 * dstar;
  }else {
    h0 = 1.5 * Minb(t,npoly+1);
  }

  if ( is.nan(h0) ){
    if ( kernel == "gauss" ){
      h0 = 0.2 * r;
    }else{
      stop("The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n")
    }
  }
  h0 <- min(h0,r)
  q = (r/(4*h0))^(1/9);
  bwCandidates = sort(q^(0:9)*h0) ;

  idx =  pracma::uniq(t)$n; # pracma

  k0_candidates <- list('quar' = 0.9375,  'epan' = 0.7500, 'rect' = 0.5000,
                        'gausvar' = 0.498677, 'gausvar1' = 0.598413,  'gausvar2' = 0.298415, 'other' = 0.398942)
  if( any(names(k0_candidates) == kernel)){
    k0 = as.numeric(k0_candidates[kernel])
  } else {
    k0 =  as.numeric(k0_candidates$other)
  }
  gcvScores <- c()
  # Get the corresponding GCV scores
  for(i in 1:length(bwCandidates)){
    newmu = fdapace::Lwls1D(bwCandidates[i], kernel_type=kernel, npoly=npoly, nder=nder, xin = t,yin= y, win = rep(1,length(y)),xout= sort(unique(t)))[idx]
    cvsum = sum((newmu -y)^2 )
    gcvScores[i] =cvsum/(1-(r*k0)/(N*bwCandidates[i]))^2
  }

  # If no bandwith gives a finite gcvScore increase the candidate bandwith and retry on a finer grid
  if(all((is.infinite(gcvScores)))){
    bwCandidates = seq( max(bwCandidates), r, length.out = 2*length(bwCandidates))
    for(i in 1:length(bwCandidates)){
      newmu = fdapace::Lwls1D(bwCandidates[i], kernel_type =kernel, npoly=npoly, nder=nder, xin = t,yin= y, win = rep(1,length(y)), xout= sort(unique(t)))[idx]
      cvsum = sum((newmu -y)^2 )
      gcvScores[i] =cvsum/(1-(r*k0)/(N*bwCandidates[i]))^2
    }
  }

  # If the problem persist we clearly have too sparse data
  if(all((is.infinite(gcvScores)))){
    stop("The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n")
  }

  bInd = which(gcvScores == min(gcvScores));
  bScr = gcvScores[bInd][1]
  bOpt = max(bwCandidates[bInd]);

  if( bOpt == r){
    warning("data is too sparse, optimal bandwidth includes all the data!You may want to change to Gaussian kernel!\n")
  }
  bOptList <- list( 'bOpt' = bOpt, 'bScore' = bScr)
  return( bOptList)
}


# This function is used to find the minimum bandwidth choice
# where the local window contains at least "numPoints" points
# Input x  : n x 1 vector
# Input numPoints: an integer specifying the number of points in a local window
# for local weighted constant, numPoints is at least 1
# for local weighted linear, numPoints is at least 2
# Output b: the minimum bandwidth choice for vector x
Minb <- function(x, numPoints){

  n = length(x);
  if( (numPoints<1) || (numPoints > n) ){
    warning("Invalid number of minimum points specified\n")
    return(NaN)
  }

  gridPts <- sort(unique(x))
  distNN1 <- max(diff(gridPts, lag=numPoints))

  return(distNN1)
}

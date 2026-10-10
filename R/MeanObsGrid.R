# Mean function estimates at all observation grids
MeanObsGrid <- function(Lt, Ly, kern, bw){
  xinraw <- unlist(Lt)
  yinraw <- unlist(Ly)

  #time in ascending order
  xin <- sort(xinraw)
  yin <- yinraw[order(xinraw)]

  xout <- unique(xin)
  meanest <- fdapace::Lwls1D(bw, kern, xin = xin, yin = yin, xout = xout)

  ret <- list()
  ret$ObsGrid <- xout
  ret$mean <- meanest

  return(ret)
}

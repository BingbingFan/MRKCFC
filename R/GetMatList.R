# Computing the list of matrices via Nadaraya-Watson estimator
GetMatList <- function(Lt, Ly, kern, bw){
  # n: the sample size
  n <- length(Lt)

  # Xlist: the list of n matrices. The i-th matrix contains the observations of the i-th subject
  # and the Nadaraya-Watson estimators of the other subjects at the observation time for
  # the i-th subject.
  Xlist <- list()
  for(i in 1:n){
    ti <- Lt[[i]]
    mi <- length(ti)
    Xlist[[i]] <- matrix(Ly[[i]], ncol = mi)
    for(j in ((1:n)[-i])){
      ks <- kader::nadwat(x = Lt[[i]], dataX = Lt[[j]], dataY = Ly[[j]], K = kernfun(kern), h = bw)
      Xlist[[i]] <- rbind(Xlist[[i]], ks)
    }
  }
  return(Xlist)
}

# Norm of vector x
normq <- function(x, q){
  if(q != 0){
    sum(abs(x)^q)^(1/q)
  }
  else{
    length(which(x != 0))
  }
}

# x is a vector that can be obtained from regularly discretizing some function f.
# This function standardize x to ensure that the integral of f is 1.
unitize <- function(x, grid = 1){
  L2 <- normq(x, 2)^(2)
  x/sqrt(L2 * grid)
}

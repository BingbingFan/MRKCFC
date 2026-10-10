#' @title Kernel Functions
#'
#' @description Some common-used kernel functions.
#'
#' @param type A \code{character} denoting the kernel type; 'epan'(Epanechnikov), 'unif'(Uniform), 'quar'(Quartic), 'gauss'(Gaussian).
#'
#' @return The corresponding kernel function.
#' @export
#'
#' @examples
#' x <- seq(-2, 2, 0.01)
#' par(mfrow = c(2,2))
#' plot(x, kernfun("epan")(x), type = "l", main = "Epanechnikov")
#' plot(x, kernfun("unif")(x), type = "l", main = "Uniform")
#' plot(x, kernfun("quar")(x), type = "l", main = "Quartic")
#' plot(x, kernfun("gauss")(x), type = "l", main = "Gaussian")
#' par(mfrow = c(1,1))
kernfun <- function(type){
  if(type == "epan"){
    return(function(x){(abs(x) < 1) * 3/4 * (1 - x^2)})
  }
  if(type == "unif"){
    return(function(x){(abs(x) < 1) * 1/2})
  }
  if(type == "quar"){
    return(function(x){(abs(x) < 1) * 15/16 * (1 - x^2)^2})
  }
  if(type == "gauss"){
    return(function(x){exp(-x^2/2)/sqrt(2*pi)})
  }
  stop("Kernel function was not set correctly!")
}

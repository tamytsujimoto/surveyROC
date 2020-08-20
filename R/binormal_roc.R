#' Binormal ROC curve
#'
#' Compute the ROC curve based on  binomal data.
#'
#' @param x numeric vector between 0 and 1
#' @param mean0 numeric value of mean for class label 0 (healthy)
#' @param mean1 numeric value of mean for class label 1 (diseased)
#' @param sigma0 numeric value of standard deviation for class label 0 (healthy)
#' @param sigma1 numeric value of standard deviation for class label 1 (diseased)
#'
#' @return a numeric vector with same dimension as x with the binormal ROC values for specified values x
#'
#' @export


binormal_roc <- function(x, mean0, mean1, sigma0, sigma1){

  a = (mean1-mean0)/sigma1
  b = sigma0/sigma1

  return(pnorm(a + b*qnorm(x)))

}

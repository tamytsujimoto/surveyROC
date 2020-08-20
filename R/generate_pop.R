#' Generate population from mixture of normal distribution
#'
#' Generate population containing healthy and diseased subpopulation generated from normal distributions
#'
#' @param N total size of population
#' @param mean0 numeric value of mean for class label 0 (healthy)
#' @param mean1 numeric value of mean for class label 1 (diseased)
#' @param sigma0 numeric value of standard deviation for class label 0 (healthy)
#' @param sigma1 numeric value of standard deviation for class label 1 (diseased)
#' @param p1 proportion of diseased observations in the population
#'
#' @return a dataframe containing
#' \itemize{
#'   \item X: continuous marker values
#'   \item D: class labels (0: healthy; 1: diseased)
#' }
#'
#' @export


generate_pop <- function(N, mean0, mean1, sigma0, sigma1, p1) {

  # Generate diseased indicator, diseased measurement and non-diseased measurement
  D <- rbinom(N, 1, p1)
  X0 <- rnorm(N, mean0, sigma0)
  X1 <- rnorm(N, mean1, sigma1)

  # Compute measurement
  X <- (1-D)*X0 +D*X1

  return(data.frame(X = X, D = D))
}

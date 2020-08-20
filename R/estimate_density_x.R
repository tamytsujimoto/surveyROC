#' Function to estimate survey weighted density
#'
#' Computes the bandwidth and estimate Gaussian kernel density estimation
#'
#' @param x numeric value for the point of interest
#' @param y numeric vector of observations from the distribution whose density is to be estimated
#' @param weight numeric vector of sampling weights
#'
#' @return a dataframe containing
#' \itemize{
#'   \item f_hat: Estimated density at point x specified
#'   \item bandwidth: Computed bandwidth
#' }
#'
#' @export

estimate_density_x <- function(x, # point of interest
                               y, # vector containing data
                               weight = NULL # vector of weights
){

  if(!is.null(weight)){

    # Computing bandwith
    N_hat <- sum(weight)
    n <- length(y)
    b <- .79*(quantile(y, .75) - quantile(y, .25))*n^(-1/5)

    f_hat <- sum(weight*dnorm((x-y)/b)/b)/N_hat
  } else{

    # Computing bandwith
    N_hat <- length(y)
    n <- length(y)
    b <- .79*(quantile(y, .75) - quantile(y, .25))*n^(-1/5)

    f_hat <- sum(dnorm((x-y)/b)/b)/N_hat
  }


  return(data.frame(f_hat = f_hat,
                    bandwidth = b))

}

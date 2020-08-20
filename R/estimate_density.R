#' Function to estimate survey weighted density
#'
#' Computes the bandwidth and estimate Gaussian kernel density estimation for a vector of desired points
#'
#' @param x vector of points of interest for density estimation
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

estimate_density <- function(x,
                             y,
                             weight = NULL){

  # Estimating density
  dens <-
    x %>%
    map_dfr(estimate_density_x,
            y = y,
            weight = weight) %>%
    mutate(x = x) %>%
    select(x, f_hat, bandwidth)

  return(dens)

}

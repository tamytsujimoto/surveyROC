#' Estimate mu_pi1
#'
#' Function to estimate mu_pi1 using Hájek estimator
#'
#' @param weight vector of sampling weights
#'
#' @return numerical value for the mu_pi1
#'
#' @export

estimate_mu_pi1 <- function(weight){

  N_hat <- sum(weight)
  mu_pi1_hat <- sum(weight^2)/N_hat - 1

  return(mu_pi1_hat)

}

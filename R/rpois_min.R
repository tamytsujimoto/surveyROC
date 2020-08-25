#' Random generator of Poisson distribution with mininum value
#'
#' @param min minimum allowed for cluster size
#' @param lambda numeric value for expected cluster size for poisson distr
#'
#' @export


rpois_min <- function(min, # minimum allowed for cluster size
                      lambda # expected cluster size for poisson distr
){
  size <- 0
  while(size < min) size <- rpois(1, lambda)
  return(size)

}

#' Generate cluster size with specified mininum size allowed
#'
#' @param n_cluster numeric value for number of clusters to be generated
#' @param min minimum allowed for cluster size
#' @param lambda numeric value for expected cluster size for poisson distr
#'
#' @export


generate_clust_size <- function(n_cluster,
                                min,
                                lambda){

  size <- replicate(n_cluster, rpois_min(min, lambda))
  return(size)

}

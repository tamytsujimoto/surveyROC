#' Generate stratified population from mixture of normal distribution
#'
#' Generate stratified population containing healthy and diseased subpopulation generated from normal distributions
#'
#' @param N total size of population
#' @param freq_strata numeric vector containing frequency of strata
#' @param mean0_strat numeric vector of mean for class label 0 (healthy) for all strata
#' @param mean1_strat numeric vector of mean for class label 1 (diseased) for all strata
#' @param sigma0_strat numeric vector of standard deviation for class label 0 (healthy) for all strata
#' @param sigma1_strat numeric vector of standard deviation for class label 1 (diseased) for all strata
#' @param p1 proportion of diseased observations in the population
#'
#' @return Dataframe containing
#' \itemize{
#'   \item strata: strata identifier
#'   \item X: continuous marker values
#'   \item D: class labels (0: healthy; 1: diseased)
#' }
#'
#' @examples
#' generate_pop_strat(N = 1000,
#'                    freq_strata = rep(.25, 4),
#'                    mean0_strat = rep(0, 4),
#'                    mean1_strat = 1:4,
#'                    sigma0_strat = rep(1, 4),
#'                    sigma1_strat = rep(1, 4),
#'                    p1 = .25)
#'
#' @export

generate_pop_strat <- function(N,
                               freq_strata,
                               mean0_strat,
                               mean1_strat,
                               sigma0_strat,
                               sigma1_strat,
                               p1) {

  N_strata <- length(freq_strata)
  size_strata <- ceiling(N*freq_strata)

  pop <-
    purrr::pmap(list(as.list(size_strata),
                     as.list(mean0_strat),
                     as.list(mean1_strat),
                     as.list(sigma0_strat),
                     as.list(sigma1_strat),
                     replicate(N_strata, p1, simplify = FALSE)),
                .f = generate_pop) %>%
    dplyr::bind_rows(.id = 'strata')

  return(pop)
}



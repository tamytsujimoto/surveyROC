#' Binormal ROC curve for stratified data
#'
#' Compute the combined binormal ROC curve for stratified data
#'
#' @param x numeric vector between 0 and 1
#' @param freq_strata numeric vector containing frequency of strata
#' @param mean0_strat numeric vector of mean for class label 0 (healthy) for all strata
#' @param mean1_strat numeric vector of mean for class label 1 (diseased) for all strata
#' @param sigma0_strat numeric vector of standard deviation for class label 0 (healthy) for all strata
#' @param sigma1_strat numeric vector of standard deviation for class label 1 (diseased) for all strata
#' @param by_strata set to TRUE if ROC curve per strata is desired
#'
#' @return a dataframe containing
#' \itemize{
#'   \item strata: if by_strata = TRUE
#'   \item x: values specified for x
#'   \item roc: values of binormal roc for each value of specified x
#' }
#'
#' @export

binormal_roc_strat <- function(x,
                               freq_strata,
                               mean0_strat,
                               mean1_strat,
                               sigma0_strat,
                               sigma1_strat,
                               by_strata = FALSE # set to TRUE if roc curve per strata is desired
){

  # Computing number of strata
  N_strata <- length(freq_strata)

  # Computing binormal ROC curve per strata
  roc_strata <-
    purrr::pmap(list(replicate(N_strata, x, simplify = FALSE),
              as.list(mean0_strat),
              as.list(mean1_strat),
              as.list(sigma0_strat),
              as.list(sigma1_strat)),
         .f = binormal_roc)

  if(by_strata == FALSE){

    # Computing theoretical ROC for population
    roc <-
      roc_strata %>%
      purrr::map2(as.list(freq_strata),
           .f = function(x,y) return(x*y)) %>%
      purrr::reduce(.f = `+`) %>%
      data.frame(x = x,
                 roc = .)

  } else{
    # Aggretaing theoretical ROC curves by strata in data.frame
    roc <-
      roc_strata %>%
      purrr::map(.f = function(y) data.frame(x = x,
                                      roc = y)) %>%
      dplyr::bind_rows(.id = 'strata') %>%
      mutate(strata = as.numeric(strata))
  }

  return(roc)

}

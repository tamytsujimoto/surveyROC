
#' Title
#'
#' @param pop
#' @param N_psu
#' @param N_sample
#' @param sample_size_psu
#' @param sample_size_ssu
#' @param grid
#'
#' @return
#' @export
#'
#' @examples

simulate_sample_tscs_strat <- function(pop,
                                       N_psu,
                                       N_sample,
                                       sample_size_psu,
                                       sample_size_ssu,
                                       grid = NULL) {


  # Generating and stacking samples
  sim <-
    replicate(N_sample,
              simulate_cum_tscs_strat2(pop,
                                      N_psu,
                                      sample_size_psu,
                                      sample_size_ssu,
                                      grid = grid), simplify = FALSE) %>%
    bind_rows(.id = 'sample')



  return(sim)


}


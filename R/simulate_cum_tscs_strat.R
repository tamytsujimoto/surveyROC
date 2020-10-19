#' Title
#'
#' @param pop
#' @param N_psu
#' @param sample_size_psu
#' @param sample_size_ssu
#' @param grid
#'
#' @return
#' @export
#'
#' @examples

simulate_cum_tscs_strat <- function(pop,
                                    N_psu,
                                    sample_size_psu,
                                    sample_size_ssu,
                                    grid = NULL){


  # Evaluating the number of necessary observations to be generated
  psu_min <- min(sample_size_psu)
  n_psu <- c(psu_min, diff(sample_size_psu))
  n_scenario <- length(sample_size_psu)

  # Generating and stacking samples
  sample_stack <-
    purrr::pmap_dfr(list(as.list(n_psu),
                         as.list(sample_size_ssu)),
                    .f = function(x,y) generate_tscs_strat(population = pop,
                                                           N_psu = N_psu,
                                                           n_psu = x,
                                                           n_ssu = y),
                    .id = 'sample_id') %>%
    mutate(wt1 = N_cluster*cluster_size/(sample_size_psu[1]*sample_size_ssu[1]),
           wt2 = N_cluster*cluster_size/(sample_size_psu[2]*sample_size_ssu[2]),
           wt3 = N_cluster*cluster_size/(sample_size_psu[3]*sample_size_ssu[3]))

  # Computing ROC for cumulative samples
  roc_cum <-
    as.list(seq(1,n_scenario)) %>%
    map(.f = function(x) as.character(seq(1,x))) %>%
    map(.f = select_sample, sample_stack = sample_stack) %>%
    map_dfr(estimate_survey_roc, grid = grid, var_weight = 'weight', ci = TRUE, .id = 'scenario')


  return(roc_cum)


}

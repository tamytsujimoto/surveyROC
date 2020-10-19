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

# Second attempt
simulate_cum_tscs_strat2 <- function(pop,
                                    N_psu,
                                    sample_size_psu,
                                    sample_size_ssu,
                                    grid = NULL){


  # Evaluating the number of necessary observations to be generated
  psu_max <- max(sample_size_psu)
  ssu_max <- max(sample_size_ssu)
  n_scenario <- length(sample_size_psu)

  # Generating and stacking samples
  sample_max <- generate_tscs_strat(population = pop,
                                    N_psu = N_psu,
                                    n_psu = psu_max,
                                    n_ssu = ssu_max) %>%
    select(strata, cluster_id, cluster_size, X, D)

  N_psu_max <- as.vector(table(sample_max$strata))

  sample_stack <-
    purrr::pmap_dfr(list(as.list(sample_size_psu),
                         as.list(sample_size_ssu)),
                    .f = function(x,y) generate_tscs_strat(population = sample_max,
                                                           N_psu = N_psu,
                                                           n_psu = x,
                                                           n_ssu = y),
                    .id = 'sample_id')


  # Computing ROC for cumulative samples
  roc_cum <-
    as.list(1:n_scenario) %>%
    map(.f = function(x) filter(sample_stack, sample_id %in% x)) %>%
    map_dfr(estimate_survey_roc, grid = grid, var_weight = 'weight', ci = TRUE, .id = 'scenario')


  return(roc_cum)


}

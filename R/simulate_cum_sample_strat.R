


simulate_cum_sample_strat <- function(pop,
                                      sample_size,
                                      grid,
                                      var_weight = 'weight') {

  # Generating and stacking samples
  sample_stack <-
    purrr::map(as.list(sample_size),
               .f = generate_srs_strat,
               pop_strat = pop)

  # Computing ROC for cumulative samples
   roc_cum <-
     sample_stack %>%
     map_dfr(estimate_survey_roc, grid = grid, var_weight = 'weight', ci = TRUE, .id = 'scenario')

  return(roc_cum)
}

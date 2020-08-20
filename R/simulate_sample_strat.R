#' Simulate stratified Simple Random Samples from a population
#'
#' Function to generate stratified SRS from a given population and compute ROC curve for each sample drawn
#'
#' @param pop dataframe containing population information (strata, X, D)
#' @param N_sample numeric value for number of samples to be drawn
#' @param sample_size sample size of samples to be drawn
#' @param grid numeric vector of points for roc curve from 0 to 1
#' @param var_weight string value indicating sampling weight variable name
#'
#' @return a dataframe containing containing survey weighted point and interval estimates for ROC curve for each
#' sample drawn from the specified population
#' @export




simulate_sample_strat <- function(pop,
                                  N_sample,
                                  sample_size,
                                  grid,
                                  var_weight = 'weight') {

  as.list(rep(sample_size, N_sample)) %>%
    map(generate_srs_strat, pop = pop) %>% # Generating stratified SR samples
    map_dfr(estimate_survey_roc, grid = grid, var_weight = var_weight, ci = TRUE, .id = 'sample')

}

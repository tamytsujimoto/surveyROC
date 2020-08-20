#' Simulation function for stratified SRS
#'
#' Generates stratified populations and draw stratified SRS from each population.
#' Computes survey weighted ROC curve for each sample drawn form each of the genereated populations
#'
#' @param N_pop numeric value for number of stratified population to be generated
#' @param pop_size numeric value for population size of populations generated
#' @param param_pop data frame containing population parameters (prob, mu0, mu1, sigma0, sigma1)
#' @param N_sample numeric value for number of samples to be drawn
#' @param sample_size numeric value for sample size of samples to be drawn
#' @param grid numeric vector of points for roc curve from 0 to 1
#'
#' @return a dataframe containing containing survey weighted point and interval estimates for ROC curve for each
#' sample drawn from each of the generated populations
#'
#' @export

simulate_strat <- function(N_pop,
                           pop_size,
                           param_pop,
                           prop_disease,
                           N_sample,
                           sample_size,
                           grid = seq(0, 1, by=0.1)){

  # Generating N_pop populations
  pop <-
    replicate(N_pop,
              generate_pop_strat(N = pop_size,
                                 freq_strata = param_pop$prob,
                                 mean0_strat = param_pop$mu0,
                                 mean1_strat = param_pop$mu1,
                                 sigma0_strat = param_pop$sigma0,
                                 sigma1_strat = param_pop$sigma1,
                                 p1 = prop_disease),
              simplify = FALSE)

  # Computing finite population ROC
  tpr_pop <-
    pop %>%
    map_dfr(estimate_survey_roc,
            grid = grid,
            .id = 'population') %>%
    select(population, fpr, tpr_pop = tpr)

  # For each population, generating N_sample samples and computing ROC
  sim <-
    pop %>%
    map_dfr(simulate_sample_strat,
            N_sample = N_sample,
            sample_size = sample_size,
            grid = grid,
            .id = 'population') %>%
    left_join(tpr_pop, by = c('fpr', 'population'))

  return(sim)

}

#' Simulation function for stratified Two Stage Cluster Sampling
#'
#' Generates stratified populations and draw stratified SRS from each population.
#' Computes survey weighted ROC curve for each sample drawn form each of the genereated populations
#'
#' @param N_pop numeric value for number of stratified population to be generated
#' @param N_psu numeric vector of number os psu's (clusters) per strata
#' @param cluster_size_strat list of cluster sizes per strata
#' @param param_pop data frame containing population parameters (prob, mu0, mu1, sigma0, sigma1)
#' @param prop_disease numeric value for the proportion of diseased in population
#' @param N_sample numeric value for number of samples to be drawn from each population
#' @param sample_size_psu numeric value for number of psu selected per strata
#' @param sample_size_ssu numeric value for number of ssu selected per selected psu
#' @param grid numeric vector of points for roc curve from 0 to 1
#'
#' @return a dataframe containing containing survey weighted point and interval estimates for ROC curve for each
#' sample drawn from each of the generated populations
#'
#' @export

simulate_tscs_strat <- function(N_pop,
                                N_psu,
                                cluster_size_strat,
                                param_pop,
                                prop_disease,
                                N_sample,
                                sample_size_psu,
                                sample_size_ssu,
                                grid = seq(0, 1, by=0.1)){



  # Generating N_pop populations
  pop <-
    replicate(N_pop,
              generate_pop_clust_strat(cluster_size = cluster_size_strat,
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
    map_dfr(simulate_sample_tscs_strat,
            N_psu = N_psu,
            N_sample = N_sample,
            sample_size_psu = sample_size_psu,
            sample_size_ssu = sample_size_ssu,
            grid = grid,
            .id = 'population') %>%
    left_join(tpr_pop, by = c('fpr', 'population'))

  return(sim)

}


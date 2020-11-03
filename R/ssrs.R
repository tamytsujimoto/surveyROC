#' Generate population from mixture of normal distribution
#'
#' Generate population containing healthy and diseased subpopulation generated from normal distributions
#'
#' @param N total size of population
#' @param mean0 numeric value of mean for class label 0 (healthy)
#' @param mean1 numeric value of mean for class label 1 (diseased)
#' @param sigma0 numeric value of standard deviation for class label 0 (healthy)
#' @param sigma1 numeric value of standard deviation for class label 1 (diseased)
#' @param p1 proportion of diseased observations in the population
#'
#' @return a dataframe containing
#' \itemize{
#'   \item X: continuous marker values
#'   \item D: class labels (0: healthy; 1: diseased)
#' }
#'
#' @export

generate_pop <- function(N, mean0, mean1, sigma0, sigma1, p1) {

  # Generate diseased indicator, diseased measurement and non-diseased measurement
  D <- rbinom(N, 1, p1)
  X0 <- rnorm(N, mean0, sigma0)
  X1 <- rnorm(N, mean1, sigma1)

  # Compute measurement
  X <- (1-D)*X0 +D*X1

  return(data.frame(X = X, D = D))
}

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
                     replicate(N_strata, p1, simplify = FALSE)), # List of parameters for each strata
                .f = generate_pop) %>% # Generating population for each strata
    dplyr::bind_rows(.id = 'strata') # Combining all results in df

  return(pop)
}

#' Generate Simple Random Sample (SRS) from population
#'
#' Generate SRS from population and compute the sampling weight
#'
#' @param pop population from where sample will be drawn
#' @param sample_size desired sample size
#'
#' @return a dataframe containing
#' \itemize{
#'   \item X: continuous marker values
#'   \item D: class labels (0: healthy; 1: diseased)
#'   \item weight: SRS weight
#'   \item pop_size: population size
#' }
#'
#' @export

generate_srs <- function(pop, sample_size){

  pop_size <- dim(pop)[1] # Population size
  weight <- pop_size/sample_size # Computing sampling weight for SRS

  sample <-
    pop %>%
    slice_sample(n = sample_size) %>% # sampling sample_size observations
    mutate(weight = weight,
           pop_size = pop_size)


  return(sample)

}

#' Generate Stratified Simple Random Sample from population
#'
#' Generate stratified SRS and compute the sampling weight
#'
#' @param pop_strat stratified population from where sample will be drawn
#' @param sample_size desired final sample size
#'
#' @return a dataframe containing
#' \itemize{
#'   \item strata: strata identifier
#'   \item strata_size: size of strata in the population
#'   \item weight: Stratified SRS weight
#'   \item X: continuous marker values
#'   \item D: class labels (0: healthy; 1: diseased)
#' }
#'
#' @export

generate_srs_strat <- function(pop_strat, sample_size){

  # Determining number of strata
  N <- dim(pop_strat)[1]
  N_strata <- length(table(pop_strat$strata))
  sample_size_strata <- ceiling(sample_size/N_strata)

  pop_nested <-
    pop_strat %>%
    group_by(strata) %>%
    nest()

  sample <-
    pop_nested$data %>%
    map(.f = generate_srs,
        sample_size = sample_size_strata) %>%
    bind_rows(.id = 'strata') %>%
    mutate(freq_strata = pop_size/N) %>%
    select(strata, strata_size = pop_size, freq_strata, weight, X, D)

  return(sample)
}

#' Simulate stratified Simple Random Samples from a population
#'
#' Function to generate stratified SRS from a given population and compute survey weighted
#' ROC curve for each sample drawn
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

  # Generating and stacking samples
  sim <-
    replicate(N_sample,
              simulate_cum_sample_strat(pop = pop,
                                        sample_size = sample_size,
                                        grid = grid,
                                        var_weight = var_weight),
              simplify = FALSE) %>%
    bind_rows(.id = 'sample')

  return(sim)

}

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


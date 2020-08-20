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

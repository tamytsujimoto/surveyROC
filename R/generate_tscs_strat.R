#' Generate a Stratified Two-Stage Cluster sampling
#'
#' @param population dataframe containing population
#' @param N_psu numeric value for number of psu's per strata in the population
#' @param n_psu numeric value for number of psu to be drawn from each strata
#' @param n_ssu numeric value for number of ssu to be drawn from each psu
#'
#' @export

generate_tscs_strat <- function(population,
                                N_psu,
                                n_psu,
                                n_ssu){

  N_strata = length(table(population$strata))
  N_psu_df = data.frame(strata = 1:N_strata,
                        N_cluster = N_psu)

  # Sampling psu
  sample_psu <-
    N_psu %>%
    map(.f = function(x) seq(1:x)) %>% # creating index per strata
    map(.f = function(x) data.frame(cluster_id = sample(x, size = n_psu))) %>% # sampling n_psu per strata
    dplyr::bind_rows(.id = 'strata') %>% # aggregating in one dataframe
    mutate(strata = as.numeric(strata),
           cluster_id = as.numeric(cluster_id)) %>%
    left_join(population, by = c('strata', 'cluster_id')) # bringing information from selected clusters

  # Sampling SSU
  final_sample <-
    sample_psu %>%
    group_by(strata, cluster_id) %>%
    slice_sample(n = n_ssu) %>%
    ungroup %>%
    left_join(N_psu_df, by = 'strata') %>%
    mutate(weight = ifelse(cluster_size < n_ssu, N_cluster/n_psu, N_cluster*cluster_size/(n_psu*n_ssu))) # computing weights

  return(final_sample)


}

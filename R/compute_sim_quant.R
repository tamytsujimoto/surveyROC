#' Simulation Results
#'
#' Function to compute relative bias, empirical standard deviation and coverage probability
#' referent to finite population and super population parameters
#'
#' @param sim simulation object
#' @param param_pop data frame containing finite population parameters (prob, mu0, mu1, sigma0, sigma1)
#'
#' @return a dataframe containing containing:
#' \itemize{
#'   \item fpr: False positive rate (equals to grid, if specified)
#'   \item tpr_th: True positive rate
#'   \item rel_bias_pop: cutoffs for continuous marker
#'   \item rel_bias_th: cutoffs for continuous marker
#'   \item emp_sd_roc: cutoffs for continuous marker
#'   \item rel_bias_pop: cutoffs for continuous marker
#'   \item rel_bias_pop: cutoffs for continuous marker
#' }
#' @export

compute_sim_quant <- function(sim, param_pop) {

  fpr <- unique(sim$fpr)

  # Computing theoretical ROC curve for population
  tpr_th <-
    binormal_roc_strat(fpr,
                       freq_strata = param_pop$prob,
                       mean0_strat = param_pop$mu0,
                       mean1_strat = param_pop$mu1,
                       sigma0_strat = param_pop$sigma0,
                       sigma1_strat = param_pop$sigma1) %>%
    select(fpr = x, tpr_th = roc)

  # Merging theoretical ROC curve and computing quantities of interest
  sim %>%
    left_join(tpr_th, by = 'fpr') %>%
    mutate(rel_bias_pop = (tpr - tpr_pop)/tpr_pop,
           rel_bias_th = (tpr - tpr_th)/tpr_th,
           ci_cov_pop = ifelse(tpr_pop >= tp_ci_l & tpr_pop <= tp_ci_u, 1, 0),
           ci_cov_th = ifelse(tpr_th >= tp_ci_l & tpr_th <= tp_ci_u, 1, 0)) %>%
    group_by(fpr) %>%
    summarise(tpr_th = unique(tpr_th),
              rel_bias_pop = mean(rel_bias_pop)*100,
              rel_bias_th = mean(rel_bias_th)*100,
              emp_sd_roc = sd(tpr),
              asy_sd_roc = mean(tpr_sd),
              ci_cov_pop = mean(ci_cov_pop)*100,
              ci_cov_th = mean(ci_cov_th)*100)

}

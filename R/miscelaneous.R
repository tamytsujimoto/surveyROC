#' Binormal ROC curve
#'
#' Compute the ROC curve based on  binormal model
#'
#' @param x numeric vector between 0 and 1
#' @param mean0 numeric value of mean for class label 0 (healthy)
#' @param mean1 numeric value of mean for class label 1 (diseased)
#' @param sigma0 numeric value of standard deviation for class label 0 (healthy)
#' @param sigma1 numeric value of standard deviation for class label 1 (diseased)
#'
#' @return a numeric vector with same dimension as x with the binormal ROC values for specified values x
#'
#' @export


binormal_roc <- function(x, mean0, mean1, sigma0, sigma1){

  a = (mean1-mean0)/sigma1
  b = sigma0/sigma1

  return(pnorm(a + b*qnorm(x)))

}

#' Binormal ROC curve for stratified data
#'
#' Compute the combined binormal ROC curve for stratified data
#'
#' @param x numeric vector between 0 and 1
#' @param freq_strata numeric vector containing frequency of strata
#' @param mean0_strat numeric vector of mean for class label 0 (healthy) for all strata
#' @param mean1_strat numeric vector of mean for class label 1 (diseased) for all strata
#' @param sigma0_strat numeric vector of standard deviation for class label 0 (healthy) for all strata
#' @param sigma1_strat numeric vector of standard deviation for class label 1 (diseased) for all strata
#' @param by_strata set to TRUE if ROC curve per strata is desired
#'
#' @return a dataframe containing
#' \itemize{
#'   \item strata: if by_strata = TRUE
#'   \item x: values specified for x
#'   \item roc: values of binormal roc for each value of specified x
#' }
#'
#' @export

binormal_roc_strat <- function(x,
                               freq_strata,
                               mean0_strat,
                               mean1_strat,
                               sigma0_strat,
                               sigma1_strat,
                               by_strata = FALSE # set to TRUE if roc curve per strata is desired
){

  # Computing number of strata
  N_strata <- length(freq_strata)

  # Computing binormal ROC curve per strata
  roc_strata <-
    purrr::pmap(list(replicate(N_strata, x, simplify = FALSE),
                     as.list(mean0_strat),
                     as.list(mean1_strat),
                     as.list(sigma0_strat),
                     as.list(sigma1_strat)),
                .f = binormal_roc)

  if(by_strata == FALSE){

    # Computing theoretical ROC for population
    roc <-
      roc_strata %>%
      purrr::map2(as.list(freq_strata),
                  .f = function(x,y) return(x*y)) %>%
      purrr::reduce(.f = `+`) %>%
      data.frame(x = x,
                 roc = .)

  } else{
    # Aggretaing theoretical ROC curves by strata in data.frame
    roc <-
      roc_strata %>%
      purrr::map(.f = function(y) data.frame(x = x,
                                             roc = y)) %>%
      dplyr::bind_rows(.id = 'strata') %>%
      mutate(strata = as.numeric(strata))
  }

  return(roc)

}

#' Binormal Area Under Curve (AUC) for ROC curve
#'
#' Compute the AUC for ROC curve based on  binormal model
#'
#' @param mean0 numeric value of mean for class label 0 (healthy)
#' @param mean1 numeric value of mean for class label 1 (diseased)
#' @param sigma0 numeric value of standard deviation for class label 0 (healthy)
#' @param sigma1 numeric value of standard deviation for class label 1 (diseased)
#'
#' @return a numeric vector with same dimension as x with the binormal ROC values for specified values x
#'
#' @export


binormal_auc <- function(mean0, mean1, sigma0, sigma1){

  a <- (mean1-mean0)/sigma1
  b <- sigma0/sigma1
  auc <- pnorm(a/sqrt(1+b^2))

  return(auc)

}


#' Binormal likelihood ratio
#'
#' Compute the AUC for ROC curve based on  binormal model
#'
#' @param c numeric vector cutoff in real line
#' @param mean0 numeric value of mean for class label 0 (healthy)
#' @param mean1 numeric value of mean for class label 1 (diseased)
#' @param sigma0 numeric value of standard deviation for class label 0 (healthy)
#' @param sigma1 numeric value of standard deviation for class label 1 (diseased)
#'
#' @return a numeric vector with same dimension as x with the likelihood ratio for binormal ROC model
#'
#' @export


binormal_lr <- function(c, mean0, mean1, sigma0, sigma1){

  a <- (mean1-mean0)/sigma1
  b <- sigma0/sigma1

  lr <- b*exp(-0.5*(c-mean1)^2/sigma1^2 + 0.5*(c-mean0)^2/sigma0^2)

  return(lr)

}

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
           ci_cov_pop = ifelse(tpr_pop >= tp_ci_l_fin & tpr_pop <= tp_ci_u_fin, 1, 0),
           ci_cov_th = ifelse(tpr_th >= tp_ci_l_sup & tpr_th <= tp_ci_u_sup, 1, 0)) %>%
    group_by(scenario, fpr) %>%
    summarise(tpr_th = unique(tpr_th),
              rel_bias_pop = mean(rel_bias_pop)*100,
              rel_bias_th = mean(rel_bias_th)*100,
              emp_sd_roc = sd(tpr),
              asy_sd_roc_sup = mean(tpr_sd_sup),
              asy_sd_roc_fin = mean(tpr_sd_fin),
              ci_cov_pop = mean(ci_cov_pop)*100,
              ci_cov_th = mean(ci_cov_th)*100) %>%
    pivot_wider(id_cols = c(scenario, fpr, tpr_th), names_from = scenario,
                values_from = -c(scenario, fpr, tpr_th))

}


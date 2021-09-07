#' Function to estimate survey weighted density
#'
#' Computes the bandwidth and estimate Gaussian kernel density estimation
#'
#' @param x numeric value for the point of interest
#' @param y numeric vector of observations from the distribution whose density is to be estimated
#' @param weight numeric vector of sampling weights
#'
#' @return a dataframe containing
#' \itemize{
#'   \item f_hat: Estimated density at point x specified
#'   \item bandwidth: Computed bandwidth
#' }
#'
#' @export

estimate_density_x <- function(x, # point of interest
                               y, # vector containing data
                               weight = NULL # vector of weights
){

  if(!is.null(weight)){

    # Computing bandwith
    N_hat <- sum(weight)
    n <- length(y)
    range <- Hmisc::wtd.quantile(y, weights = weight, probs = .75) - Hmisc::wtd.quantile(y, weights = weight, probs = .25)

    b <- .79*(range)*N_hat^(-1/5)

    f_hat <- sum(weight*dnorm((x-y)/b)/b)/N_hat
  } else{

    # Computing bandwith
    N_hat <- length(y)
    n <- length(y)
    b <- .79*(quantile(y, .75) - quantile(y, .25))*n^(-1/5)

    f_hat <- sum(dnorm((x-y)/b)/b)/N_hat
  }


  return(data.frame(f_hat = f_hat,
                    bandwidth = b))

}


#' Function to estimate survey weighted density
#'
#' Computes the bandwidth and estimate Gaussian kernel density estimation for a vector of desired points
#'
#' @param x vector of points of interest for density estimation
#' @param y numeric vector of observations from the distribution whose density is to be estimated
#' @param weight numeric vector of sampling weights
#'
#' @return a dataframe containing
#' \itemize{
#'   \item f_hat: Estimated density at point x specified
#'   \item bandwidth: Computed bandwidth
#' }
#'
#' @export

estimate_density <- function(x,
                             y,
                             weight = NULL){

  # Estimating density
  dens <-
    x %>%
    map_dfr(estimate_density_x,
            y = y,
            weight = weight) %>%
    mutate(x = x) %>%
    select(x, f_hat, bandwidth)

  return(dens)

}

#' Estimate mu_pi1
#'
#' Function to estimate mu_pi1 using Hájek estimator
#'
#' @param weight vector of sampling weights
#'
#' @return numerical value for the mu_pi1
#'
#' @export

estimate_mu_pi1 <- function(weight){

  N_hat <- sum(weight)
  mu_pi1_hat <- sum(weight^2)/N_hat - 1

  return(mu_pi1_hat)

}


#' Estimate ROC curve for complex survey data
#'
#' Function to estimate ROC curve (FPR, TPR) for complex survey data
#'
#' @param data data frame containing continuous marker X and class label D
#' @param grid optional numeric vector of points for roc curve from 0 to 1
#' @param var_weight string variable name for sampling weight. If specified, should be included in data
#' @param ci Logical, if true, will calculate asymptotic confidence interval for the TPF
#'
#' @return a dataframe containing
#' \itemize{
#'   \item fpr: False positive rate (equals to grid, if specified)
#'   \item tpr: True positive rate
#'   \item cutoff: cutoffs for continuous marker
#' }
#'
#' @export

fD1 <- function(t) {
  f1 <- dnorm(t, param_pop$mu1[1], param_pop$sigma1[1])
  f2 <- dnorm(t, param_pop$mu1[2], param_pop$sigma1[2])
  f3 <- dnorm(t, param_pop$mu1[3], param_pop$sigma1[3])
  f4 <- dnorm(t, param_pop$mu1[4], param_pop$sigma1[4])
  f5 <- dnorm(t, param_pop$mu1[5], param_pop$sigma1[5])

  param_pop$prob[1]*f1 +
    param_pop$prob[2]*f2 +
    param_pop$prob[3]*f3 +
    param_pop$prob[4]*f4 +
    param_pop$prob[5]*f5
}

roc_var_svy <- function(roc_data,
                        design,
                        grid = NULL){

  if(is.null(grid)) grid <- seq(0,1,.01)

  roc_svy <-
    roc_data %>%
    mutate(grid = cut(fpr, c(-Inf,grid), include.lowest = TRUE, labels = grid)) %>%
    group_by(grid) %>%
    summarise(fpr = max(fpr),
              tpr = max(tpr),
              cutoff = min(cutoff))

  cutoff <- roc_svy$cutoff
  fp <- roc_svy$fpr
  tp <- roc_svy$tpr

  n <- dim(design$variables)[1]
  wt <- 1/(design$prob)
  mu_pi1 <- estimate_mu_pi1(wt)
  p1 <- svymean(~D, design)[1]
  sampling_frac <- n/sum(wt)

  dsgnD1 <- subset(design, D == 1)
  dsgnD0 <- subset(design, D == 0)

  f0 =
    svysmooth(~X, design = dsgnD0)$X %>%
    as_tibble() %>%
    mutate(aux = cut(x, c(-Inf,cutoff), include.lowest = TRUE)) %>%
    group_by(aux) %>%
    filter(row_number() == n()) %>%
    select(aux, f0 = y)

  f1 =
    svysmooth(~X, design = dsgnD1)$X %>%
    as_tibble() %>%
    mutate(aux = cut(x, c(-Inf,cutoff), include.lowest = TRUE)) %>%
    group_by(aux) %>%
    filter(row_number() == n()) %>%
    select(aux, f1 = y)

  ratio <-
    f1 %>%
    full_join(f0, by = 'aux') %>%
    mutate(ratio = f1/f0)

  roc_svy <-
    roc_svy %>%
    mutate(aux = cut(cutoff, c(-Inf,cutoff), include.lowest = TRUE)) %>%
    left_join(ratio, by = "aux") %>%
    mutate(#ratio_true = fD1(cutoff)/dnorm(cutoff),
      vsvy_sup = (1/n)*(sampling_frac*(1+mu_pi1))*((1-p1)^(-1)*ratio^2*fpr*(1-fpr) + p1^(-1)*tpr*(1-tpr)),
      vsvy_fin = (1/n)*(sampling_frac*(mu_pi1))*((1-p1)^(-1)*ratio^2*fpr*(1-fpr) + p1^(-1)*tpr*(1-tpr))) %>%
    select(-aux)

  return(roc_svy)

}

estimate_survey_roc <- function(design,
                                grid = NULL,
                                var_weight = NULL,
                                ci = FALSE) {

  data <- design$variables
  D <- data$D
  X <- data$X
  X.order <- order(X, decreasing = TRUE)

  n <- dim(data)[1]
  wt <- rep(1, n)
  if(!is.null(var_weight)) wt <- data[[var_weight]]

  TTT <- X[X.order]
  TPF <- cumsum(wt[X.order]*(D[X.order] == 1))/sum(wt*(D == 1))
  FPF <- cumsum(wt[X.order]*(D[X.order] == 0))/sum(wt*(D == 0))

  dups <- rev(duplicated(rev(TTT)))
  tp <- TPF[!dups]
  fp <- FPF[!dups]
  cutoffs <- TTT[!dups]

  tp <- ifelse(tp > 1, 1, tp)
  fp <- ifelse(fp > 1, 1, fp)

  df_roc <-
    data.frame(fpr = c(0,fp),
                       tpr = c(0,tp),
                       cutoff = c(Inf, cutoffs))

  # CONFIDENCE INTERVAL

  if(ci) {

    df_roc <-
      roc_var_svy(df_roc, design, grid) %>%
      mutate(ci_ll_sup = tpr - 1.96*sqrt(vsvy_sup),
             ci_ul_sup = tpr + 1.96*sqrt(vsvy_sup),
             ci_ll_fin = tpr - 1.96*sqrt(vsvy_fin),
             ci_ul_fin = tpr + 1.96*sqrt(vsvy_fin))

  }

  if(ci == FALSE & !is.null(grid)) {

    df_roc <-
      df_roc %>%
      mutate(grid = cut(fpr, c(-Inf,grid), include.lowest = TRUE, labels = grid)) %>%
      group_by(grid) %>%
      summarise(fpr = max(fpr),
                tpr = max(tpr),
                cutoff = min(cutoff))

  }

  return(df_roc)
}


#' Funtion to compute and plot the survey roc
#'
#' @param data
#' @param var_X
#' @param var_D
#' @param var_weight
#' @param ci
#'
#' @return
#' @export
#'
#' @examples

plot_survey_roc <- function(data, var_X, var_D, var_weight = NULL, ci = FALSE){

  roc_data <-
    data %>%
    select_if(names(.) %in% c(var_X, var_D, var_weight)) %>%
    rename(X = .data[[var_X]],
           D = .data[[var_D]]) %>%
    filter(complete.cases(.))

  # COMPUTING ROC

  roc <-
    roc_data %>%
    estimate_survey_roc(var_weight = var_weight)

  # COMPUTING ROC CI

  if(ci){
    roc_ci <-
      roc_data %>%
      estimate_survey_roc(var_weight = var_weight,
                          grid = seq(0,1,.01),
                          ci = TRUE)
  }

  # PLOTTING CURVE

  p <-
    if(ci){
      roc_ci %>%
        ggplot(aes(x = fpr, y = tpr)) +
        #geom_step(aes(x = fpr, y = tp_ci_l)) +
        #geom_step(aes(x = fpr, y = tp_ci_u)) +
        geom_ribbon(aes(ymin = tp_ci_l, ymax = tp_ci_u, fill = "95% Confidence Interval"), alpha = .5) +
        geom_step(aes(x = fpr, y = tpr, colour = "Survey weighted"), data = roc) +
        theme_bw() +
        labs(x = 'FPR', y = 'TPR')+
        theme(legend.position = "top") +
        scale_colour_manual("", values=c("Unweighted" = "red", "Survey weighted" = "black")) +
        scale_fill_manual("", values=c("95% Confidence Interval" = "grey"))
    } else{
      roc %>%
        ggplot(aes(x = fpr, y = tpr)) +
        geom_step() +
        theme_bw()
    }


  result <-
    if(ci) list(roc, roc_ci, p)
  else list(roc, p)

  return(result)

}

#' Function to compute Area Under the Curve using trapezoid method
#'
#' @param fpr
#' @param tpr
#'
#' @return
#' @export
#'
#' @examples

compute_auc <- function(fpr, tpr){

  auc <- 0
  for (i in 2:length(fpr)) {
    auc <- auc + 0.5 * (fpr[i] - fpr[i - 1]) * (tpr[i] + tpr[i - 1])
  }

  return(auc)

}







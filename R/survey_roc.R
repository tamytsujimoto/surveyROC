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
    b <- .79*(quantile(y, .75) - quantile(y, .25))*n^(-1/5)

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


estimate_survey_roc <- function(data,
                                grid = NULL,
                                var_weight = NULL,
                                ci = FALSE) {

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

  # CONFIDENCE INTERVAL

  if(!is.null(grid)){

    aux <- cut(fp, grid, include.lowest = TRUE)
    tp <- as.vector(tapply(tp, aux, max))
    fp <- grid[-1]
    cutoffs <- as.vector(tapply(cutoffs, aux, min))
  }

  df <- data.frame(fpr = c(0,fp),
                   tpr = c(0,tp),
                   cutoff = c(Inf, cutoffs))

  if(!is.null(var_weight) & ci){

    # Estimating densities
    dens0 <- estimate_density(x = cutoffs, y = X[D == 0], weight = wt[D == 0])
    dens1 <- estimate_density(x = cutoffs, y = X[D == 1], weight = wt[D == 1])
    f0 <- dens0$f_hat
    f1 <- dens1$f_hat

    # Estimating sampling fraction, mu_pi1 and p1
    p1 <- sum(D*wt)/sum(wt)
    sampling_frac <- n/sum(wt)
    mu_pi1 <- estimate_mu_pi1(wt)

    # Computing asymptotic standard deviation and 95% CI superpop
    tp_var_super = (sampling_frac*(1+mu_pi1))*((1-p1)^(-1)*(f1/f0)^2*fp*(1-fp) + p1^(-1)*tp*(1-tp))
    tp_var_finite = (sampling_frac*(mu_pi1))*((1-p1)^(-1)*(f1/f0)^2*fp*(1-fp) + p1^(-1)*tp*(1-tp))

    tp_sd_super = sqrt(tp_var_super/n)
    tp_sd_finite = sqrt(tp_var_finite/n)

    tp_ci_l_super = tp - qnorm(.975)*tp_sd_super
    tp_ci_u_super = tp + qnorm(.975)*tp_sd_super
    tp_ci_l_finite = tp - qnorm(.975)*tp_sd_finite
    tp_ci_u_finite = tp + qnorm(.975)*tp_sd_finite

    #
    df <- cbind(df,
                data.frame(
                  f0 = c(NA, f0),
                  bandwidth0 = c(NA, dens0$bandwidth),
                  f1 = c(NA, f1),
                  bandwidth1 = c(NA, dens1$bandwidth),
                  tpr_sd_sup = c(0, tp_sd_super),
                  tpr_sd_fin = c(0, tp_sd_finite),
                  tp_ci_l_sup = c(0, tp_ci_l_super),
                  tp_ci_u_sup = c(0, tp_ci_u_super),
                  tp_ci_l_fin = c(0, tp_ci_l_finite),
                  tp_ci_u_fin = c(0, tp_ci_u_finite)
                ))
  }

  return(df)
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







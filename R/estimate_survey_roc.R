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

    # Computing asymptotic standard deviation and 95% CI
    tp_var = (sampling_frac*(1+mu_pi1))*((1-p1)^(-1)*(f1/f0)^2*fp*(1-fp) + p1^(-1)*tp*(1-tp))
    tp_sd = sqrt(tp_var/n)
    tp_ci_l = tp - qnorm(.975)*tp_sd
    tp_ci_u = tp + qnorm(.975)*tp_sd

    #
    df <- cbind(df,
                data.frame(
                  f0 = c(NA, f0),
                  bandwidth0 = c(NA, dens0$bandwidth),
                  f1 = c(NA, f1),
                  bandwidth1 = c(NA, dens1$bandwidth),
                  tpr_sd = c(0, tp_sd),
                  tp_ci_l = c(0, tp_ci_l),
                  tp_ci_u = c(0, tp_ci_u)))
  }

  return(df)
}



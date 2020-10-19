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




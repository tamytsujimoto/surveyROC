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




#' Title
#'
#' @param id
#' @param sample_stack
#'
#' @return
#' @export
#'
#' @examples

select_sample <- function(id, sample_stack){

  sample <-
    sample_stack %>%
    filter(sample_id %in% id) %>%
    select(strata, cluster_id, cluster_size, X, D, N_cluster, weight)

  return(sample)
}

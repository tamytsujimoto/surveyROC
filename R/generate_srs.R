#' Generate Simple Random Sample (SRS) from population
#'
#' Generate SRS from population and compute the sampling weight
#'
#' @param pop population from where sample will be drawn
#' @param sample_size desired sample size
#'
#' @return a dataframe containing
#' \itemize{
#'   \item X: continuous marker values
#'   \item D: class labels (0: healthy; 1: diseased)
#'   \item weight: SRS weight
#'   \item pop_size: population size
#' }
#'
#' @export

generate_srs <- function(pop, sample_size){

  pop_size <- dim(pop)[1]
  weight <- pop_size/sample_size

  sample <-
    pop %>%
    slice_sample(n = sample_size) %>% # sampling sample_size observations
    mutate(weight = weight,
           pop_size = pop_size)


  return(sample)

}










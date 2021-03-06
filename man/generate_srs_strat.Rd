% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssrs.R
\name{generate_srs_strat}
\alias{generate_srs_strat}
\title{Generate Stratified Simple Random Sample from population}
\usage{
generate_srs_strat(pop_strat, sample_size)
}
\arguments{
\item{pop_strat}{stratified population from where sample will be drawn}

\item{sample_size}{desired final sample size}
}
\value{
a dataframe containing
\itemize{
  \item strata: strata identifier
  \item strata_size: size of strata in the population
  \item weight: Stratified SRS weight
  \item X: continuous marker values
  \item D: class labels (0: healthy; 1: diseased)
}
}
\description{
Generate stratified SRS and compute the sampling weight
}

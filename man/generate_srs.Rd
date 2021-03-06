% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssrs.R
\name{generate_srs}
\alias{generate_srs}
\title{Generate Simple Random Sample (SRS) from population}
\usage{
generate_srs(pop, sample_size)
}
\arguments{
\item{pop}{population from where sample will be drawn}

\item{sample_size}{desired sample size}
}
\value{
a dataframe containing
\itemize{
  \item X: continuous marker values
  \item D: class labels (0: healthy; 1: diseased)
  \item weight: SRS weight
  \item pop_size: population size
}
}
\description{
Generate SRS from population and compute the sampling weight
}

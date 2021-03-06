% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssrs.R
\name{generate_pop}
\alias{generate_pop}
\title{Generate population from mixture of normal distribution}
\usage{
generate_pop(N, mean0, mean1, sigma0, sigma1, p1)
}
\arguments{
\item{N}{total size of population}

\item{mean0}{numeric value of mean for class label 0 (healthy)}

\item{mean1}{numeric value of mean for class label 1 (diseased)}

\item{sigma0}{numeric value of standard deviation for class label 0 (healthy)}

\item{sigma1}{numeric value of standard deviation for class label 1 (diseased)}

\item{p1}{proportion of diseased observations in the population}
}
\value{
a dataframe containing
\itemize{
  \item X: continuous marker values
  \item D: class labels (0: healthy; 1: diseased)
}
}
\description{
Generate population containing healthy and diseased subpopulation generated from normal distributions
}

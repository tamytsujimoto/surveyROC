% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survey_roc.R
\name{estimate_density}
\alias{estimate_density}
\title{Function to estimate survey weighted density}
\usage{
estimate_density(x, y, weight = NULL)
}
\arguments{
\item{x}{vector of points of interest for density estimation}

\item{y}{numeric vector of observations from the distribution whose density is to be estimated}

\item{weight}{numeric vector of sampling weights}
}
\value{
a dataframe containing
\itemize{
  \item f_hat: Estimated density at point x specified
  \item bandwidth: Computed bandwidth
}
}
\description{
Computes the bandwidth and estimate Gaussian kernel density estimation for a vector of desired points
}

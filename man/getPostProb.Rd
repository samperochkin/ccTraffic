% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_functions.R
\name{getPostProb}
\alias{getPostProb}
\title{Compute posterior probabilities from a posterior distribution (fitted model).}
\usage{
getPostProb(fit, beta_names, quantiles, M = 10000)
}
\arguments{
\item{fit}{A fitted case-crossover model (obtained via `fitModel(model)`).}

\item{beta_names}{Names of beta parameters concerned.}

\item{quantiles}{Quantiles at which to compute the probabilities.}

\item{M}{Number of replicates to draw for the posterior (to approximate the probability, the more the better).}
}
\value{
A vector of (approximate) probabilities.
}
\description{
Compute posterior probabilities from a posterior distribution (fitted model).
}

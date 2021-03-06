% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_adjust_batch.R
\name{diagnostic_adjust_batch}
\alias{diagnostic_adjust_batch}
\title{Diagnostic visualization for adj_batch function}
\usage{
diagnostic_adjust_batch(feature_abd, feature_abd_adj, var_batch, gamma_hat,
  gamma_star, output)
}
\arguments{
\item{feature_abd}{original feature-by-sample matrix of abundances
(proportions or counts).}

\item{feature_abd_adj}{feature-by-sample matrix of batch-adjusted feature
abundances, with covariate effects retained and scales consistent with
original abundance matrix.}

\item{var_batch}{the batch variable (should be a factor).}

\item{gamma_hat}{estimated per feature-batch gamma parameters.}

\item{gamma_star}{shrinked per feature-batch gamma parameters}

\item{output}{output file name}
}
\value{
(invisbly) the ggplot2 plot object
}
\description{
Diagnostic visualization for adj_batch function
}
\keyword{internal}

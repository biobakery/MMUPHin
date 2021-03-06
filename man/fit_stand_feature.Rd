% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_adjust_batch.R
\name{fit_stand_feature}
\alias{fit_stand_feature}
\title{Fit lm and standardize all features}
\usage{
fit_stand_feature(s_data, design, l_ind)
}
\arguments{
\item{s_data}{feature-by-sample matrix of abundances (proportions or
counts).}

\item{design}{design matrix.}

\item{l_ind}{list of indicator matrices, as returned by construct_ind.}
}
\value{
list of two componet: the standardized feature abundance matrix, and
a list of per-feature standardization fits.
}
\description{
Fit lm and standardize all features
}
\keyword{internal}

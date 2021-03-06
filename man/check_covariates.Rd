% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_lm_meta.R
\name{check_covariates}
\alias{check_covariates}
\title{Check covariates}
\usage{
check_covariates(data_covariates, batch)
}
\arguments{
\item{data_covariates}{data frame of covariates.}

\item{batch}{batch variable.}
}
\value{
vector of indicators per batch for if/which covariates can be fitted 
within the batches
}
\description{
Check covariates
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{check_feature_abd}
\alias{check_feature_abd}
\title{Check feature abundance table}
\usage{
check_feature_abd(feature_abd)
}
\arguments{
\item{feature_abd}{feature-by-sample matrix of abundances (proportions or
counts).}
}
\value{
returns an error if any of the check fails. Otherwise either "counts"
or "proportions"
}
\description{
Given a feature abundance table, make sure that a) it has no missing values,
b) all values are non-negative, c) it is either proportions (all no greater
than 1) or counts (all integers).
}
\keyword{internal}

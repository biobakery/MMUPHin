% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{normalize_features}
\alias{normalize_features}
\title{Normalize feature abundance table (modified from Maaslin2)}
\usage{
normalize_features(features, normalization = "NONE", pseudo_count = 0)
}
\arguments{
\item{features}{feature-by-sample matrix of abundances (proportions or
counts).}

\item{normalization}{normalization method.}

\item{pseudo_count}{pseudo count to be added to feature_abd.}
}
\value{
normalized abundance table.
}
\description{
Normalize feature abundance table (modified from Maaslin2)
}
\keyword{internal}

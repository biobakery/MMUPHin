% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_adjust_batch.R
\name{it_sol}
\alias{it_sol}
\title{Iteratively solve for one feature's shrinked location and scale parameters}
\usage{
it_sol(s_data, g_hat, d_hat, g_bar, t2, a, b, control)
}
\arguments{
\item{s_data}{the feature's standardized abundances.}

\item{g_hat}{the feature's location parameter frequentist estimations.}

\item{d_hat}{the feature's scale parameter frequentist estimations.}

\item{g_bar}{EB estimation of location hyper parameters.}

\item{t2}{EB estimation of location hyper parameters.}

\item{a}{EB estimation of scale hyper parameters.}

\item{b}{EB estimation of scale hyper parameters.}

\item{control}{list of control parameters}
}
\value{
matrix of shrinked location and scale parameters.
}
\description{
Iteratively solve for one feature's shrinked location and scale parameters
}
\keyword{internal}

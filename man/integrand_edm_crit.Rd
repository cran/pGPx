% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integrand_edm_crit.R
\name{integrand_edm_crit}
\alias{integrand_edm_crit}
\title{Integrand of the distance in measure criterion}
\usage{
integrand_edm_crit(
  x,
  E,
  model,
  Thresh,
  batchsize,
  alpha,
  predE,
  predx = NULL,
  precalc.data = NULL
)
}
\arguments{
\item{x}{vector of dimension \eqn{d} representing the \eqn{ith} point where to compute the criterion}

\item{E}{matrix of dimension \eqn{d*(i-1)} containing the previously optimized simulation points}

\item{model}{km model}

\item{Thresh}{threshold selected for excursion set}

\item{batchsize}{number of simulation points}

\item{alpha}{value of Vorob'ev threshold}

\item{predE}{list containing the posterior mean and standard deviation at E}

\item{predx}{list containing the posterior mean and standard deviation at x}

\item{precalc.data}{list result of \link[KrigInv]{precomputeUpdateData} with \code{model} and \code{x}.}
}
\value{
the value of the integrand at \eqn{x}
}
\description{
Computes the integrand of the distance in measure criterion.
}
\references{
Azzimonti D. F., Bect J., Chevalier C. and Ginsbourger D. (2016). Quantifying uncertainties on excursion sets under a Gaussian random field prior. SIAM/ASA Journal on Uncertainty Quantification, 4(1):850–874.

Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
}

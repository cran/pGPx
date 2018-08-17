#' @title Distance in measure criterion
#' @description Computes the distance in measure criterion. To be used in optimization routines.
# Input:
#' @param x vector of dimension \eqn{d} representing the point where to compute the criterion
#' @param other.points Vector giving the other batchsize-1 points at which one wants to evaluate the criterion
#' @param integration.points p*d matrix of points for numerical integration in the X space.
#' @param integration.weights Vector of size p corresponding to the weights of these integration points.
#' @param intpoints.oldmean Vector of size p corresponding to the kriging mean at the integration points.
#' @param intpoints.oldsd Vector of size p corresponding to the kriging standard deviation at the integration points.
#' @param precalc.data list result of \link[KrigInv]{precomputeUpdateData} with \code{model} and \code{x}.
#' @param model km model
#' @param threshold threshold selected for excursion set
#' @param batchsize number of simulation points
#' @param alpha value of Vorob'ev threshold
#' @param current.crit Current value of the criterion
#' @return the value of the expected distance in measure criterion at \eqn{x},\code{other.points}.
#' @references Azzimonti D. F., Bect J., Chevalier C. and Ginsbourger D. (2016). Quantifying uncertainties on excursion sets under a Gaussian random field prior. SIAM/ASA Journal on Uncertainty Quantification, 4(1):850–874.
#'
#' Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
#' @export
edm_crit2 <- function (x, other.points, integration.points, integration.weights = NULL,
                      intpoints.oldmean, intpoints.oldsd, precalc.data, model,
                      threshold, batchsize, alpha, current.crit)
{
  x.complete <- c(x,other.points)
  return(edm_crit(
    x = x.complete, integration.points = integration.points, integration.weights = integration.weights,
    intpoints.oldmean = intpoints.oldmean,intpoints.oldsd = intpoints.oldsd,precalc.data = precalc.data,
    model = model,threshold = threshold,batchsize=batchsize,alpha=alpha,
    current.crit=current.crit)
  )
}

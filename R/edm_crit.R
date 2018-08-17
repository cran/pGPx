#' @title Distance in measure criterion
#' @description Computes the distance in measure criterion.
# Input:
#' @param x vector of dimension \eqn{d} representing the \eqn{ith} point where to compute the criterion
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
#' @return the value of the expected distance in measure criterion at \eqn{x}
#' @references Azzimonti D. F., Bect J., Chevalier C. and Ginsbourger D. (2016). Quantifying uncertainties on excursion sets under a Gaussian random field prior. SIAM/ASA Journal on Uncertainty Quantification, 4(1):850–874.
#'
#' Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
#' @export
edm_crit <- function (x, integration.points, integration.weights = NULL,
          intpoints.oldmean, intpoints.oldsd, precalc.data, model,
          threshold, batchsize, alpha, current.crit)
{
  d <- model@d
  n <- model@n
  X.new <- matrix(x, nrow = d)
  mindist <- Inf
  tp1 <- c(as.numeric(t(model@X)), x)
  for (i in 1:batchsize) {
    xx <- X.new[, i]
    tp2 <- matrix(tp1 - as.numeric(xx), ncol = d, byrow = TRUE)^2
    mysums <- sqrt(rowSums(tp2))
    mysums[n + i] <- Inf
    mindist <- min(mindist, mysums)
  }

  if (!identical(colnames(integration.points), colnames(model@X)))
    colnames(integration.points) <- colnames(model@X)





  if ((mindist > 1e-05) ) {
    X.new <- t(X.new)
    predE <-KrigInv::predict_nobias_km(object = model, newdata = X.new,
                                       type = "UK", se.compute = TRUE, cov.compute = TRUE)

    predx <-list(mean= intpoints.oldmean, sd= intpoints.oldsd)

    result<-integrand_edm_crit(x=t(integration.points),E=X.new,model=model,
                               Thresh = threshold,batchsize = batchsize,
                               alpha = alpha,predE = predE,predx = predx,precalc.data = precalc.data)

    if (is.null(integration.weights)) {
      crit <- mean(result)
    }
    else crit <- sum(result * integration.weights)
  }
  else crit <- current.crit + 0.01
  return(crit)
}

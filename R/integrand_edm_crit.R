#' @title Integrand of the distance in measure criterion
#' @description Computes the integrand of the distance in measure criterion.
# Input:
#' @param x vector of dimension \eqn{d} representing the \eqn{ith} point where to compute the criterion
#' @param E matrix of dimension \eqn{d*(i-1)} containing the previously optimized simulation points
#' @param model km model
#' @param Thresh threshold selected for excursion set
#' @param batchsize number of simulation points
#' @param alpha value of Vorob'ev threshold
#' @param predE list containing the posterior mean and standard deviation at E
#' @param predx list containing the posterior mean and standard deviation at x
#' @param precalc.data list result of \link[KrigInv]{precomputeUpdateData} with \code{model} and \code{x}.
#' @return the value of the integrand at \eqn{x}
#' @references Azzimonti D. F., Bect J., Chevalier C. and Ginsbourger D. (2016). Quantifying uncertainties on excursion sets under a Gaussian random field prior. SIAM/ASA Journal on Uncertainty Quantification, 4(1):850–874.
#'
#' Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
#' @export
integrand_edm_crit = function (x, E, model, Thresh, batchsize, alpha,predE,
                               predx = NULL,precalc.data=NULL)
{

  x<-t(x)

  if(is.null(predx)){
    predx <- KrigInv::predict_nobias_km(object=model,newdata=x,type="UK",se.compute=TRUE)
  }
  evalpoint.oldmean <- predx$mean ;
  evalpoint.oldsd <- ifelse(predx$sd<.Machine$double.eps,0,predx$sd);


  if(is.null(precalc.data))
    precalc.data <- KrigInv::precomputeUpdateData(model,x)

  d <- model@d

  X.new <- matrix(t(E), nrow = d)

  if (!identical(colnames(x), colnames(model@X)))
    colnames(x) <- colnames(model@X)

  X.new <- t(X.new)
  mk <- predE$mean
  sk <- predE$sd
  newXvar <- sk * sk
  F.newdata <- predE$F.newdata
  c.newdata <- predE$c
  Sigma.r <- predE$cov
  kn = KrigInv::computeQuickKrigcov(model, x, X.new,
                           precalc.data, F.newdata, c.newdata)
  krig2 <- KrigInv::predict_update_km_parallel(newXmean = mk, newXvar = newXvar,
                                      newXvalue = mk, Sigma.r = Sigma.r, newdata.oldmean = evalpoint.oldmean,
                                      newdata.oldsd = evalpoint.oldsd, kn = kn)
  if (!is.null(krig2$error)){
    current.vorob <- pnorm((evalpoint.oldmean-Thresh)/evalpoint.oldsd)
    return(current.vorob)
  }


  gamma_n <- rowSums(krig2$lambda*kn)
  gamma_n <- ifelse(gamma_n <.Machine$double.eps,0, gamma_n)

  arg1 <- as.numeric((evalpoint.oldmean - Thresh)/evalpoint.oldsd)
  arg2 <- as.numeric((Thresh-evalpoint.oldmean)/sqrt(gamma_n))
  arg3 <- -sqrt(gamma_n)/evalpoint.oldsd

  arg1[arg1 == -Inf] <- -1000
  arg1[arg1 == Inf] <- 1000
  arg1[is.nan(arg1)] <- 1000


  arg2[arg2 == -Inf] <- -1000
  arg2[arg2 == Inf] <- 1000
  arg2[is.nan(arg2)] <- -1000*sign(arg1)


  arg3[arg3 <= -1] <- -1
  arg3[arg3 >= 1] <- 1
  arg3[is.nan(arg3)] <- 0

  term1 <- pbivnorm::pbivnorm(arg1, arg2, arg3)
  term2 <- pbivnorm::pbivnorm(-arg1, -arg2, arg3)
  crit <- term1 + term2

  return(crit)
}

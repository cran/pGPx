#' @title Simulate and interpolate
#' @description Generates \code{nsims} approximate posterior field realizations
#' at \code{interpolatepoints}. The approximate realizations are computed by
#' simulating the field only at \code{simupoints} simulation points.
# Input:
#' @param object km object
#' @param nsim numbero of simulations
#' @param simupoints simulations points, locations where the field was simulated
#' @param interpolatepoints points where posterior simulations are approximated
#' @param nugget.sim nugget to be added to simulations for numerical stability
#' @param type type of kriging model used for approximation (default Universal Kriging)
#' @return A matrix \code{nsim*interpolatepoints} containing the approximate realizations.
#' @references Azzimonti D. F., Bect J., Chevalier C. and Ginsbourger D. (2016). Quantifying uncertainties on excursion sets under a Gaussian random field prior. SIAM/ASA Journal on Uncertainty Quantification, 4(1):850–874.
#'
#' Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
#' @examples
#' ### Simulate and interpolate for a 2d example
#' if (!requireNamespace("DiceKriging", quietly = TRUE)) {
#' stop("DiceKriging needed for this example to work. Please install it.",
#'      call. = FALSE)
#' }
#' if (!requireNamespace("DiceDesign", quietly = TRUE)) {
#' stop("DiceDesign needed for this example to work. Please install it.",
#'      call. = FALSE)
#' }
#' # Define the function
#' g=function(x){
#'   return(-DiceKriging::branin(x))
#' }
#' d=2
#' # Fit OK km model
#' design<-DiceDesign::maximinESE_LHS(design = DiceDesign::lhsDesign(n=50,
#'                                                                   dimension = 2,
#'                                                                   seed=42)$design)$design
#' colnames(design)<-c("x1","x2")
#' observations<-apply(X = design,MARGIN = 1,FUN = g)
#' kmModel<-DiceKriging::km(formula = ~1,design = design,response = observations,
#'                          covtype = "matern3_2",control=list(trace=FALSE))
#' # Get simulation points
#' # Here they are not optimized, you can use optim_dist_measure to find optimized points
#' simu_points <- DiceDesign::maximinSA_LHS(DiceDesign::lhsDesign(n=100,
#'                                                                dimension = d,
#'                                                                seed=1)$design)$design
#' # obtain nsims posterior realization at simu_points
#' nsims <- 1
#' nn_data<-expand.grid(seq(0,1,,50),seq(0,1,,50))
#' nn_data<-data.frame(nn_data)
#' colnames(nn_data)<-colnames(kmModel@X)
#' approx.simu <- simulate_and_interpolate(object=kmModel, nsim = 1, simupoints = simu_points,
#'                                         interpolatepoints = as.matrix(nn_data),
#'                                         nugget.sim = 0, type = "UK")
#' \donttest{
#' ## Plot the approximate process realization
#' image(matrix(approx.simu[1,],ncol=50),
#'       col=grey.colors(20))
#' contour(matrix(approx.simu[1,],ncol=50),
#'         nlevels = 20,add=TRUE)
#'
#' }
#' @export
simulate_and_interpolate <- function(object, nsim = 1,
                                     simupoints = NULL, interpolatepoints = NULL,
                                     nugget.sim = 0,  type="UK"){

  some.simu <- simulate(object=object,nsim=nsim,newdata=simupoints,nugget.sim=nugget.sim,
                           cond=TRUE,checkNames = FALSE) # simulations of GRF with mean and cov given by predict(object=object,newdata=simupoints,type=type)

  obj <- krig_weight_GPsimu(object=object,simu_points=simupoints,krig_points=interpolatepoints)

  krig.mean.init <- matrix(obj$krig.mean.init,ncol=1)
  weights <- t(obj$Lambda.end)

  result <- matrix(nrow=nsim,ncol=nrow(interpolatepoints))

  for(i in 1:nsim){
    some.simu.i <- matrix(some.simu[i,],nrow=1)
    result[i,] <- krig.mean.init + tcrossprod(weights,some.simu.i)
  }

  return(result)

}

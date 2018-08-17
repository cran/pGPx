#' @title Choose simulation points
#' @description Selects \code{batchsize} locations where to simulate the field by minimizing the distance in measure criterion or by maximizing the integrand of the distance in measure criterion. Currently it is only a wrapper for the functions \code{max_distance_measure} and \code{max_integrand_edm}.
# Input:
#' @param model a km model
#' @param threshold threshold value
#' @param lower a \eqn{d} dimensional vector containing the lower bounds for the optimization
#' @param upper a \eqn{d} dimensional vector containing the upper bounds for the optimization
#' @param batchsize number of simulations points to find
#' @param algorithm type of algorithm used to obtain simulation points: \itemize{
#'        \item \code{"A"} minimize the full integral criterion;
#'        \item \code{"B"} maximize the integrand of the criterion.
#' }
#' @param alpha value of Vorob'ev threshold
#' @param verb an integer to choose the level of verbosity
#' @param optimcontrol a list containing optional parameters for the optimization, see \link[KrigInv]{max_sur_parallel} for more details.
#' @param integration.param a list containing parameters for the integration of the criterion A, see \link[KrigInv]{max_sur_parallel} for more details.
#' @return A list containing \itemize{
#'         \item \code{par} a matrix \code{batchsize*d} containing the optimal points
#'         \item \code{value} a vector of length \code{batchsize} with the values of the criterion at each step
#'}
#' @references Azzimonti D. F., Bect J., Chevalier C. and Ginsbourger D. (2016). Quantifying uncertainties on excursion sets under a Gaussian random field prior. SIAM/ASA Journal on Uncertainty Quantification, 4(1):850–874.
#'
#' Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
#' @examples
#' ### Compute optimal simulation points in a 2d example
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
#' design<-DiceDesign::maximinESE_LHS(design = DiceDesign::lhsDesign(n=20,
#'                                                                   dimension = 2,
#'                                                                   seed=42)$design)$design
#' colnames(design)<-c("x1","x2")
#' observations<-apply(X = design,MARGIN = 1,FUN = g)
#' kmModel<-DiceKriging::km(formula = ~1,design = design,response = observations,
#'                          covtype = "matern3_2",control=list(trace=FALSE))
#'
#' # Run optim_dist_measure, algorithm B to obtain one simulation point
#' # NOTE: the approximating process resulting from 1 simulation point
#' # is very rough and it should not be used, see below for a more principled example.
#' simu_pointsB <- optim_dist_measure(model=kmModel,threshold = -10,
#'                                   lower = c(0,0),upper = c(1,1),
#'                                   batchsize = 1,algorithm = "B")
#'
#' \dontrun{
#' # Get 75 simulation points with algorithm A
#' batchsize <- 50
#' simu_pointsA <- optim_dist_measure(model=kmModel,threshold = -10,
#'                                   lower = c(0,0),upper = c(1,1),
#'                                   batchsize = batchsize,algorithm = "A")
#'
#' # Get 75 simulation points with algorithm B
#' batchsize <- 75
#' simu_pointsB <- optim_dist_measure(model=kmModel,threshold = -10,
#'                                   lower = c(0,0),upper = c(1,1),
#'                                   batchsize = batchsize,algorithm = "B")
#' # plot the criterion value
#' critValA <-c(simu_pointsA$value,rep(NA,25))
#' par(mar = c(5,5,2,5))
#' plot(1:batchsize,critValA,type='l',main="Criterion value",ylab="Algorithm A",xlab="batchsize")
#' par(new=T)
#' plot(1:batchsize,simu_pointsB$value, axes=F, xlab=NA, ylab=NA,col=2,lty=2,type='l')
#' axis(side = 4)
#' mtext(side = 4, line = 3, 'Algorithm B')
#' legend("topright",c("Algorithm A","Algorithm B"),lty=c(1,2),col=c(1,2))
#' par(mar= c(5, 4, 4, 2) + 0.1)
#'
#' # obtain nsims posterior realization at simu_points
#' nsims <- 1
#' nn_data<-expand.grid(seq(0,1,,50),seq(0,1,,50))
#' nn_data<-data.frame(nn_data)
#' colnames(nn_data)<-colnames(kmModel@X)
#' approx.simu <- simulate_and_interpolate(object=kmModel, nsim = 1, simupoints = simu_pointsA$par,
#'                                         interpolatepoints = as.matrix(nn_data),
#'                                         nugget.sim = 0, type = "UK")
#'
#' ## Plot the approximate process realization
#' image(matrix(approx.simu[1,],ncol=50),
#'       col=grey.colors(20))
#' contour(matrix(approx.simu[1,],ncol=50),
#'         nlevels = 20,add=TRUE)
#' points(simu_pointsA$par,pch=17)
#' points(simu_pointsB$par,pch=1,col=2)
#'
#' }
#' @export
optim_dist_measure <- function(model, threshold, lower, upper, batchsize, algorithm= "B", alpha=0.5,verb=1, optimcontrol=NULL, integration.param=NULL){

  # Fix dimension
  d <- model@d


  # select criterion to optimize
  if(algorithm=="A"){

    if(is.null(optimcontrol))
      optimcontrol <- list(method="genoud",pop.size=100,print.level=verb)

    if(is.null(integration.param)){
      integcontrol <- list(distrib="sobol",n.points=1000)
      integration.param <- integration_design(integcontrol=integcontrol,d=d,
                                              lower=lower,upper=upper,model=model,T=threshold)
      integration.param$alpha <- 0.5
    }

    o<-max_distance_measure(lower=lower,upper = upper,optimcontrol = optimcontrol,batchsize = batchsize,integration.param = integration.param,T = threshold,model = model)
  }else{
    o<-max_integrand_edm(lower = lower,upper=upper,batchsize = batchsize,alpha=alpha,Thresh = threshold,model = model,verb = verb)
  }


  return(list(par=o$par, value=o$value))
}

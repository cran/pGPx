## Compute distance in measure for a set of simulation points
#' @name expDistMeasure
#' @title Compute expected distance in measure of approximate excursion set
#' @description Computes expected distance in measure between the excursion set of the approximated process and the true excursion set.
# Input:
#' @param simupoints a numeric array of size \code{batchsize*d} containing the simulation points.
#' @param model a km model
#' @param threshold threshold value
#' @param batchsize number of simulations points
#' @param integration.param a list containing parameters for the integration of the criterion A, see \link[KrigInv]{max_sur_parallel} for more details.
#' @return A positive value indicating the expected distance in measure.
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
#'
#' threshold <- -10
#'
#' # Obtain simulation point sampling from maximin LHS design
#' batchsize <- 50
#' set.seed(1)
#' mmLHS_simu_points <-  DiceDesign::maximinSA_LHS(DiceDesign::lhsDesign(n=batchsize,
#'                                                                      dimension = d,
#'                                                                      seed=1)$design)$design
#'
#'
#' # Compute expected distance in measure for approximation obtain from random simulation points
#' EDM_mmLHS <- rep(NA,batchsize)
#' integcontrol <- list(distrib="sobol",n.points=1000)
#' integration.param <- KrigInv::integration_design(integcontrol,d=d,
#'                                         lower=c(0,0),upper=c(1,1),
#'                                         model=kmModel,T=threshold)
#' integration.param$alpha <- 0.5
#' for(i in seq(1,batchsize)){
#' EDM_mmLHS[i]<-expDistMeasure( mmLHS_simu_points[1:i,],model = kmModel,
#'                              threshold = threshold,batchsize = i,
#'                              integration.param = integration.param  )
#' }
#'
#' plot(EDM_mmLHS,type='l',main="Expected distance in measure",xlab="batchsize")
#'
#'
#' \dontrun{
#' # Get optimized simulation points with algorithm B
#' simu_points <- optim_dist_measure(model=kmModel,threshold = threshold,
#'                                   lower = c(0,0),upper = c(1,1),
#'                                   batchsize = batchsize,algorithm = "B")
#' # plot the criterion value
#' plot(1:batchsize,simu_points$value,type='l',main="Criterion value")
#'
#' # Compute expected distance in measure for approximation obtained from optimized simulation points
#' EDM_optB <- rep(NA,batchsize)
#' for(i in seq(1,batchsize)){
#'   EDM_optB[i]<-expDistMeasure( simu_points$par[1:i,],model = kmModel,threshold = threshold,
#'                                  batchsize = i,integration.param = integration.param  )
#' }
#' plot(EDM_mmLHS,type='l',main="Expected distance in measure",
#'      xlab="batchsize",ylab="EDM",
#'      ylim=range(EDM_mmLHS,EDM_optB))
#' lines(EDM_optB,col=2,lty=2)
#' legend("topright",c("Maximin LHS","B"),lty=c(1,2),col=c(1,2))
#'
#' # Get optimized simulation points with algorithm A
#' simu_pointsA <- optim_dist_measure(model=kmModel,threshold = threshold,
#'                                    lower = c(0,0),upper = c(1,1),
#'                                    batchsize = batchsize,algorithm = "A")
#' # plot the criterion value
#' plot(1:batchsize,simu_pointsA$value,type='l',main="Criterion value")
#'
#' # Compute expected distance in measure for approximation obtained from optimized simulation points
#' EDM_optA <- rep(NA,batchsize)
#' for(i in seq(1,batchsize)){
#'   EDM_optA[i]<-expDistMeasure( simu_pointsA$par[1:i,],model = kmModel,threshold = threshold,
#'                                  batchsize = i,integration.param = integration.param  )
#' }
#' plot(EDM_mmLHS,type='l',main="Expected distance in measure",
#'      xlab="batchsize",ylab="EDM",
#'      ylim=range(EDM_mmLHS,EDM_optB,EDM_optA))
#' lines(EDM_optB,col=2,lty=2)
#' lines(EDM_optA,col=3,lty=3)
#' legend("topright",c("Maximin LHS","A","B"),lty=c(1,3,2),col=c(1,3,2))
#'
#' }
#' @export
expDistMeasure<- function(simupoints,model,threshold,batchsize,integration.param=NULL){


  if(is.null(integration.param)){
    integcontrol <- list(distrib="sobol",n.points=1000)
    integration.param <- integration_design(integcontrol,d=d,lower=c(0,0),upper=c(1,1),model=model,T=threshold)
    integration.param$alpha <- 0.5
  }

  integration.points <- as.matrix(integration.param$integration.points) ; d <- model@d
  integration.weights <- integration.param$integration.weights
  alpha <- integration.param$alpha


  #precalculates the kriging mean and variance on the integration points
  pred <- predict_nobias_km(object=model,newdata=integration.points,type="UK",se.compute=TRUE)
  intpoints.oldmean <- pred$mean ; intpoints.oldsd <- pred$sd; pn <- pnorm((intpoints.oldmean-threshold)/intpoints.oldsd)

  if(is.null(alpha)) alpha <- vorob_threshold(pn)
  pn_bigger_than_alpha <- (pn>alpha)+0
  pn_lower_than_alpha <- 1-pn_bigger_than_alpha

  if(is.null(integration.weights)) current.vorob <- mean(pn*pn_lower_than_alpha + (1-pn)*pn_bigger_than_alpha)
  if(!is.null(integration.weights)) current.vorob <- sum(integration.weights*(pn*pn_lower_than_alpha + (1-pn)*pn_bigger_than_alpha))
  precalc.data <- precomputeUpdateData(model,integration.points)


  distConfig<-edm_crit(x = t(simupoints),integration.points =integration.points,integration.weights = integration.weights,
                       intpoints.oldmean = intpoints.oldmean,intpoints.oldsd = intpoints.oldsd,
                       precalc.data = precalc.data,model = model,threshold = threshold,batchsize = batchsize,
                       alpha = alpha,current.crit = current.vorob )


  return(distConfig)
}

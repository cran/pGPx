#' @title Compute Distance Transform Variability
#' @description Compute the expected L^2 distance between the average distance transform and the set realizations.    If the input is the actual values of the gaussian process, compute also the random sets.
# Input:
#' @param rand.set a matrix of size \code{n.int.points}x\code{nsim} containing the excursion set realizations stored as long vectors. For example the excursion set obtained from the result of \code{\link{simulate_and_interpolate}}.
#' @param threshold threshold value
#' @param nsim number of simulations of the excursion set
#' @param n.int.points total length of the excursion set discretization. The size of the image is \code{sqrt(n.int.points)}.
#' @return A list containing\itemize{
#'     \item{\code{variance:}}{Value of the distance transform variability. The integral of \code{dvar} over the spatial domain.}
#'     \item{\code{dbar:}}{empirical distance average transform \eqn{ 1/N \sum_{i=1}^N d(x,\Gamma_i)}, a matrix of size \code{n.int.points} x \code{n.int.points}}
#'     \item{\code{dvar:}}{empirical variance of distance transform \eqn{ 1/N \sum_{i=1}^N (d(x,\Gamma_i) - dbar)^2}, a matrix of size \code{n.int.points} x \code{n.int.points}}
#'     \item{\code{alldt:}}{distance transforms for all realizations, a matrix of size \code{n.int.points} x \code{nsim}}
#'     \item{\code{naTot:}}{Total number of infinite distance transform values. These are returned in realizations where there is no excursion.}
#' }
#' @references Azzimonti D. F., Bect J., Chevalier C. and Ginsbourger D. (2016). Quantifying uncertainties on excursion sets under a Gaussian random field prior. SIAM/ASA Journal on Uncertainty Quantification, 4(1):850–874.
#'
#' Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
#'
#' Felzenszwalb, P. F. and Huttenlocher, D. P. (2012). Distance Transforms of Sampled Functions. Theory of Computing, 8(19):415-428.
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
#'
#' # obtain nsims posterior realization at simu_points
#' nsims <- 30
#' nn_data<-expand.grid(seq(0,1,,50),seq(0,1,,50))
#' nn_data<-data.frame(nn_data)
#' colnames(nn_data)<-colnames(kmModel@X)
#' approx.simu <- simulate_and_interpolate(object=kmModel, nsim = nsims, simupoints = simu_points,
#'                                         interpolatepoints = as.matrix(nn_data),
#'                                         nugget.sim = 0, type = "UK")
#' Dvar<- DTV(rand.set = approx.simu,threshold = -10,
#'                              nsim = nsims,n.int.points = 50^2)
#' \donttest{
#' image(matrix(Dvar$dbar,ncol=50),col=grey.colors(20),main="average distance transform")
#' image(matrix(Dvar$dvar,ncol=50),col=grey.colors(20),main="variance of distance transform")
#' points(design,pch=17)
#' }
#' @export
DTV = function(rand.set,threshold,nsim,n.int.points){


  if(!is.logical(rand.set))
    rand.set <- (rand.set > threshold)#(sign(PG-threshold)+1)/2

  if(ncol(rand.set)!=nsim & nrow(rand.set)==nsim)
    rand.set<-t(rand.set)

  #for each random set: compute the distance transform
  #then compute sum((dist_i - d_bar)^2)
  #and return that

  #the integration points are a grid
  matsize<-sqrt(n.int.points)
  alldt<-matrix(c(0),nrow=n.int.points,ncol=nsim)

  for(i in 1:nsim){
    mymat  <- matrix(rand.set[,i] , matsize, matsize )
    mydt  <- dtt_fast( mymat )
    alldt[,i]<-c(mydt)/matsize
  }

  # modification to exclude NA (infinite) distance tranforms. Those DT result from sets that are all false
  alldt[which(alldt>sqrt(2))]<-NA
  dbar <-  rowMeans(alldt,na.rm=TRUE)

  L2norm <- ( alldt - dbar )^2
  dvar <-  rowMeans(L2norm,na.rm=TRUE)
  naTot<-sum(is.na(alldt))
  res<-mean(dvar)

  return(list(variance=res,dbar=dbar,dvar=dvar,alldt=alldt,naTot=naTot))
}

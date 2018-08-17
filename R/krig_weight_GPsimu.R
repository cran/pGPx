#' @title Weights for interpolating simulations
#' @description Returns a list with the posterior mean and the kriging weights for simulations points.
# Input:
#' @param object km object.
#' @param simu_points simulations points, locations where the field was simulated.
#' @param krig_points points where the interpolation is computed.
#' @param T.mat a matrix (n+p)x(n+p) representing the Choleski factorization of the covariance matrix for the initial design and simulation points.
#' @param F.mat a matrix (n+p)x(fdim) representing the evaluation of the model matrix at the initial design and simulation points.
#' @return A list containing the posterior mean and the (ordinary) kriging weights for simulation points.
#' @references Azzimonti D. F., Bect J., Chevalier C. and Ginsbourger D. (2016). Quantifying uncertainties on excursion sets under a Gaussian random field prior. SIAM/ASA Journal on Uncertainty Quantification, 4(1):850–874.
#'
#' Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
#' @examples
#' ######################################################################
#' ### Compute the weights for approximating process on a 1d example
#' if (!requireNamespace("DiceKriging", quietly = TRUE)) {
#' stop("DiceKriging needed for this example to work. Please install it.",
#'      call. = FALSE)
#' }
#' if (!requireNamespace("DiceDesign", quietly = TRUE)) {
#' stop("DiceDesign needed for this example to work. Please install it.",
#'      call. = FALSE)
#' }
#' ## Create kriging model from GP realization
#' design<-DiceDesign::maximinESE_LHS(design = DiceDesign::lhsDesign(n=20,
#'                                                                  dimension = 1,
#'                                                                  seed=42)$design)$design
#' colnames(design)<-c("x1")
#' gp0 <- DiceKriging::km (formula = ~1, design = design,
#'                         response = rep (x = 0, times = nrow (design)),
#'                         covtype = "matern3_2", coef.trend = 0,
#'                         coef.var = 1, coef.cov = 0.2)
#' set.seed(1)
#' observations <- t (DiceKriging::simulate (object = gp0, newdata = design, cond = FALSE))
#'
#' # Fit OK km model
#' kmModel<-DiceKriging::km(formula = ~1,design = design,response = observations,
#'                          covtype = "matern3_2",control=list(trace=FALSE))
#'
#' # Get simulation points
#' # Here they are not optimized, you can use optim_dist_measure to find optimized points
#' set.seed(2)
#' simu_points <- matrix(runif(20),ncol=1)
#' # obtain nsims posterior realization at simu_points
#' nsims <- 10
#' set.seed(3)
#' some.simu <- DiceKriging::simulate(object=kmModel,nsim=nsims,newdata=simu_points,nugget.sim=1e-6,
#'                          cond=TRUE,checkNames = FALSE)
#' grid<-seq(0,1,,100)
#' nn_data<-data.frame(grid)
#' colnames(nn_data)<-colnames(kmModel@X)
#' pred_nn<-DiceKriging::predict.km(object = kmModel,newdata = nn_data,type = "UK")
#' obj <- krig_weight_GPsimu(object=kmModel,simu_points=simu_points,krig_points=grid)
#'
#' \donttest{
#' # Plot the posterior mean and some approximate process realizations
#' result <- matrix(nrow=nsims,ncol=length(grid))
#'
#' plot(nn_data$x1,pred_nn$mean,type='l')
#' for(i in 1:nsims){
#'    some.simu.i <- matrix(some.simu[i,],ncol=1)
#'    result[i,] <- obj$krig.mean.init + crossprod(obj$Lambda.end,some.simu.i)
#'    points(simu_points,some.simu.i)
#'    lines(grid,result[i,],col=3)
#' }
#' }
#' @export
krig_weight_GPsimu <- function( object, simu_points, krig_points,T.mat=NULL,F.mat=NULL){


  # p: number of points where we simulate at
  # object: a km object with n observation.
  # simu_points: the p points where the simulation is performed
  # krig_points: the points where the kriging mean is calculated

  d <- ncol(simu_points)
  p <- nrow(simu_points)
  krig_points<-matrix(krig_points,ncol=d)

  r <- nrow(krig_points)

  simu_points.mat <- matrix(simu_points,ncol=d)
  krig_points.mat <- matrix(krig_points,ncol=d)
  krig_points<-data.frame(krig_points)
  colnames(krig_points)<-colnames(object@X)

  # we want to compute kriging weights of all the p simulation points of simu_points
  # (and also the design points)
  # for the prediction at the r points krig_points
  # total number of weights: r (p+n)

  # we rely on UK weights equations (Chevalier Ph.D., eq (2.12))
  alldata <- rbind(object@X,simu_points)
  alldata.mat <- as.matrix(alldata)

  if(!is.null(F.mat)){
    F <- F.mat
  }else{
    F <- model.matrix(object=object@trend.formula, data = data.frame(alldata))
  }
  fx <- t(model.matrix(object=object@trend.formula, data = krig_points))
  if(!is.null(T.mat)){
    T<-T.mat
  }else{
    K <- covMatrix(object=object@covariance,X=alldata.mat)$C
    T <- chol(K)
  }
  M <- backsolve(T, F, upper.tri = TRUE,transpose = TRUE)
  inv_tF.Kinv.F <- chol2inv(chol(crossprod(M)))

  kx <- covMat1Mat2(object=object@covariance,X1=alldata.mat,X2=krig_points.mat)
  Tinv.kx <- backsolve(T, kx, upper.tri=TRUE,transpose = TRUE)
  Kinv.kx <- backsolve(T, Tinv.kx, upper.tri=TRUE)

  tmp <- fx - crossprod(F,Kinv.kx)
  member2 <- kx + F%*%crossprod(inv_tF.Kinv.F , tmp) #crossprod( t(F) ,crossprod(inv_tF.Kinv.F , tmp))
  Tinv.member2 <- backsolve(T, member2, upper.tri=TRUE,transpose = TRUE)
  Kinv.member2 <- backsolve(T, Tinv.member2, upper.tri=TRUE)

  lambda <- Kinv.member2

  n_plus_p <- nrow(lambda)
  n <- n_plus_p - p
  Lambda.init <- lambda[c(1:n),]
  Lambda.end <- lambda[c((n+1):n_plus_p),]

  krig.mean.init <- crossprod ( object@y ,  Lambda.init)


#   ## compute k_n(simu_points,krig_points)
#   k.iniz<-covMat1Mat2(object = object@covariance,X1 = alldata.mat,X2=krig_points.mat)
#   kSimu <- covMat1Mat2(object=object@covariance,X1 = alldata.mat,X2=simu_points.mat)
#   secondPart<-crossprod(kSimu,Kinv.kx)
#
#   fSimu <- t(model.matrix(object=object@trend.formula, data = data.frame(simu_points)))
#   Tinv.kSimu <- backsolve(t(T), kSimu, upper.tri=FALSE)
#   Kinv.kSimu <- backsolve(T, Tinv.kSimu, upper.tri=TRUE)
#   tmp2<-fSimu-crossprod(F,Kinv.kSimu)
#   innerPart<-crossprod(t(inv_tF.Kinv.F),tmp2)
#
#   thirdPart<-crossprod(tmp,innerPart)


  return(list (krig.mean.init=as.numeric(krig.mean.init),Lambda.end=Lambda.end  ))

}

#' @title Gradient of the weights for interpolating simulations
#' @description Returns a list with the gradients of the posterior mean and the gradient of the (ordinary) kriging weights for simulations points.
# Input:
#' @param object km object
#' @param simu_points simulations points, locations where the field was simulated.
#' @param krig_points one point where the interpolation is computed.
#' @param T.mat a matrix (n+p)x(n+p) representing the Choleski factorization of the covariance matrix for the initial design and simulation points.
#' @param F.mat a matrix (n+p)x(fdim) representing the evaluation of the model matrix at the initial design and simulation points.
#' @return A list containing the gradients of posterior mean and kriging weights for simulation points.
#' @references Azzimonti D. F., Bect J., Chevalier C. and Ginsbourger D. (2016). Quantifying uncertainties on excursion sets under a Gaussian random field prior. SIAM/ASA Journal on Uncertainty Quantification, 4(1):850–874.
#'
#' Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
#' @examples
#' ######################################################################
#' ### Compute the weights and gradient on a 2d example
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
#' set.seed(1)
#' simu_points <- matrix(runif(100*d),ncol=d)
#' # obtain nsims posterior realization at simu_points
#' nsims <- 1
#' set.seed(2)
#' some.simu <- DiceKriging::simulate(object=kmModel,nsim=nsims,newdata=simu_points,nugget.sim=1e-6,
#'                          cond=TRUE,checkNames = FALSE)
#' nn_data<-expand.grid(seq(0,1,,50),seq(0,1,,50))
#' nn_data<-data.frame(nn_data)
#' colnames(nn_data)<-colnames(kmModel@X)
#' obj<-krig_weight_GPsimu(object = kmModel,simu_points = simu_points,krig_points = as.matrix(nn_data))
#'
#' \donttest{
#' ## Plot the approximate process realization and the gradient vector field
#' k_scale<-5e-4
#' image(matrix(obj$krig.mean.init+crossprod(obj$Lambda.end,some.simu[1,]),ncol=50),
#'       col=grey.colors(20))
#' contour(matrix(obj$krig.mean.init+crossprod(obj$Lambda.end,some.simu[1,]),ncol=50),
#'         nlevels = 20,add=TRUE)
#'
#' for(c_ii in c(1,seq(10,2500,by = 64))){
#'    pp<-t(as.matrix(nn_data)[c_ii,])
#'    obj_deriv <- grad_kweights(object = kmModel,simu_points = simu_points,krig_points = pp)
#'    S_der<-obj_deriv$krig.mean.init + crossprod(obj_deriv$Lambda.end,some.simu[1,])
#'    points(x = pp[1],y = pp[2],pch=16)
#'    arrows(x0=pp[1],y0=pp[2],x1 = pp[1]+k_scale*S_der[1,1],y1=pp[2]+k_scale*S_der[2,1])
#' }
#' }
#' @export
grad_kweights<-function (object,  simu_points, krig_points,T.mat=NULL,F.mat=NULL){

  d <- ncol(simu_points)
  p <- nrow(simu_points)
  krig_points<-matrix(krig_points,ncol=d)

  r <- nrow(krig_points)

  simu_points.mat <- matrix(simu_points,ncol=d)
  krig_points.mat <- matrix(krig_points,ncol=d)


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
  # Derivative of
#  fx <- t(model.matrix(object=object@trend.formula, data = data.frame(krig_points)))
  beta <- object@trend.coef
  # Pre process formula to remove things that are not in derivative tables
  newFormula<-as.formula(gsub(")","",gsub("I(","",object@trend.formula,fixed=TRUE),fixed=TRUE))

  der.formula<-deriv(newFormula,namevec = colnames(object@X),function.arg = TRUE)
  trendDeriv<-function(x){
    result<-matrix(NA,ncol=ncol(x),nrow=nrow(x))
    for(i in seq(nrow(x)))
      result[i,]<-attr(do.call(der.formula,as.list(x[i,])),"gradient")
    return(result) #cbind(rep(0,nrow(x)),result))
  }
  fx <-  trendDeriv(krig_points.mat)

  if(!is.null(T.mat)){
    T <- T.mat
  }else{
    K <- covMatrix(object=object@covariance,X=alldata.mat)$C
    T <- chol(K)
  }
  M <- backsolve(T, F, upper.tri = TRUE,transpose = TRUE)
  inv_tF.Kinv.F <- chol2inv(chol(crossprod(M)))


  # Here we derive
#  kx <- covMat1Mat2(object=object@covariance,X1=alldata.mat,X2=krig_points.mat)
  c.newdata <- covMat1Mat2(object@covariance, X1 = alldata, X2 = krig_points.mat,
                           nugget.flag = object@covariance@nugget.flag)
  c.deriv.newdata<-covVector.dx(object = object@covariance,x = krig_points.mat,X = alldata,c = c.newdata) #cbind(rep(0,nrow(alldata)),covVector.dx(object = object@covariance,x = krig_points.mat,X = alldata,c = c.newdata))


  Tinv.kx <- backsolve(T, c.deriv.newdata, upper.tri=TRUE,transpose = TRUE)
  Kinv.kx <- backsolve(T, Tinv.kx, upper.tri=TRUE)

  tmp <- as.vector(fx) - crossprod(F,Kinv.kx)
  member2 <- c.deriv.newdata + F%*%crossprod(inv_tF.Kinv.F , tmp)
  Tinv.member2 <- backsolve(T, member2, upper.tri=TRUE,transpose = TRUE)
  Kinv.member2 <- backsolve(T, Tinv.member2, upper.tri=TRUE)

  lambda <- Kinv.member2

  n_plus_p <- nrow(lambda)
  n <- n_plus_p - p
  Lambda.init <- lambda[c(1:n),]
  Lambda.end <- lambda[c((n+1):n_plus_p),]

  krig.mean.init <- crossprod ( object@y ,  Lambda.init)

  return(list (krig.mean.init=as.numeric(krig.mean.init),Lambda.end=Lambda.end  ))
}

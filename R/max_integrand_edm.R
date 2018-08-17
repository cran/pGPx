#' @title Maximize the integrand distance in measure criterion
#' @description Optimizes the integrand of the distance in measure criterion.
# Input:
#' @param lower a \eqn{d} dimensional vector containing the lower bounds for the optimization
#' @param upper a \eqn{d} dimensional vector containing the upper bounds for the optimization
#' @param batchsize number of simulations points to find
#' @param alpha value of Vorob'ev threshold
#' @param Thresh threshold value
#' @param model a km model
#' @param verb an integer to choose the level of verbosity
#' @return A list containing \itemize{
#'         \item \code{par} a matrix \code{batchsize*d} containing the optimal points
#'         \item \code{value} a vector of length \code{batchsize} with the value of the criterion after each optimization
#'         \item \code{fcount} count of the number of criterion evaluations
#'}
#' @references Azzimonti D. F., Bect J., Chevalier C. and Ginsbourger D. (2016). Quantifying uncertainties on excursion sets under a Gaussian random field prior. SIAM/ASA Journal on Uncertainty Quantification, 4(1):850–874.
#'
#' Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
#' @export
max_integrand_edm <- function(lower, upper, batchsize,alpha=0.5,Thresh, model,verb=1){


  d <- model@d

  fun.optim <- integrand_edm_crit

  fcount<-0

  #Optimization with Optim
  vecMatch <- function(x, want) {
    out <- apply(x,1, function(x, want) isTRUE(all.equal(x, want)), want)
    any(out)
  }


  selectIP<-function(n,d,model,Thresh){
    Points<-matrix(runif(n*d),ncol = d)
    evalIP<-predict_nobias_km(object = model, newdata = as.data.frame(Points),type = "UK", se.compute = TRUE, cov.compute = FALSE,low.memory = TRUE)
    pn <- pnorm((evalIP$mean - Thresh)/evalIP$sd)

    sur<-pn*(1-pn)
    ss<-sort(sur,index.return=TRUE,decreasing = TRUE)
    #  ordInitPoints<-initPoints[ss$ix,]

    return(Points[ss$ix,])
  }

  #  ordInitPoints<-maximinLHS(batchsize*60,k = d)
  ordInitPoints<-matrix(sobol(n=batchsize*100*d,dim = d),ncol = d)

  evalOIP<-predict_nobias_km(object = model, newdata = as.data.frame(ordInitPoints),type = "UK", se.compute = TRUE, cov.compute = FALSE,low.memory = TRUE)
  pn <- pnorm((evalOIP$mean - Thresh)/evalOIP$sd)
  sur<-pn*(1-pn)

  #  ss<-sort(sur,index.return=TRUE,decreasing = TRUE)
  #  ordInitPoints<-ordInitPoints[ss$ix,]
  #  E<-ordInitPoints[1,]
  #  ordInitPoints<-ordInitPoints[2:(batchsize*d+1),]
  #  shuffle<-sample(batchsize*d)
  #  ordInitPoints<-ordInitPoints[shuffle,] + c(rnorm(batchsize*d,mean = 0,sd = 0.01),rnorm(batchsize*d,mean = 0,sd = 0.01))
  pointsUsed<-0
  iiMax<-which.max(sur)
  E<-ordInitPoints[iiMax,]
  values<-rep(NA,batchsize)
  indexes<-seq(batchsize*d*100)
  ss<-sample(indexes[-iiMax],size = batchsize*max(d,6),prob = sur[-iiMax])
  ordInitPoints<-matrix(ordInitPoints[ss,],ncol = d)
  #batchsize optimisations in dimension d
  for (i in 1:batchsize){

    if(verb>0)
      cat("\n ----------------\nStarting optimization of point number ",i, " of ",batchsize,"\n ----------------\n\n")

    predE <- predict_nobias_km(object = model, newdata = as.data.frame(t(matrix(t(E), nrow = d))),type = "UK", se.compute = TRUE, cov.compute = TRUE)

    o<-optim(par=ordInitPoints[pointsUsed+1,],fn=fun.optim,gr=NULL,alpha=alpha,E=E, model=model, Thresh=Thresh, batchsize=i,
             predE=predE,method="L-BFGS-B",lower=lower,upper=upper,control=list(trace=verb,lmm=15,factr=1e4,fnscale=-1))
    fcount<-fcount+o$counts[1]
    pointsUsed<-pointsUsed+1
    aa<-0
    while(o$convergence!=0){
      aa<-aa+1
      if(verb>0)
        cat("\n ---------------- \n Failed convergence, new optim attempt \n Batchsize: ",batchsize,",iteraction: ",i," (",aa,")\n ----------------\n")
      #  if(aa==1)
      oIP<- ordInitPoints[pointsUsed+1,]        # selectIP(50*d,d,model,Thresh)
      pointsUsed<-pointsUsed+1
      o<-optim(par=oIP,fn=fun.optim,gr=NULL,alpha=alpha,E=E, model=model, Thresh=Thresh, batchsize=i,
               predE=predE,method="L-BFGS-B",lower=lower,upper=upper,control=list(trace=verb,lmm=15,factr=1e4,fnscale=-1))
      fcount<-fcount+o$counts[1]
    }
    if(i>1){
      while(vecMatch(t(matrix(E,nrow = d)),o$par)){
        aa<-aa+1
        if(verb>0)
          cat("\n ---------------- \nConvergence to existing point, new optim attempt \n Batchsize: ",batchsize,",iteraction: ",i," (",aa,")\n ----------------\n")
        #    if(aa==1)
        oIP<- ordInitPoints[pointsUsed+1,]       # selectIP(50*d,d,model,Thresh)
        pointsUsed<-pointsUsed+1
        o<-optim(par=oIP,fn=fun.optim,gr=NULL,alpha=alpha,E=E, model=model, Thresh=Thresh, batchsize=i,
                 predE=predE,method="L-BFGS-B",lower=lower,upper=upper,control=list(trace=verb,lmm=15,factr=1e4,fnscale=-1))
        fcount<-fcount+o$counts[1]
      }
    }

    E <- c(E,as.numeric(o$par))
    values[i]<-o$value
    if(verb>0)
      cat("\nPoint number ",i, " of ",batchsize," optimized\n ----------------\n\n")
  }
  o$par <- t(matrix(E,nrow=d)); colnames(o$par) <- colnames(model@X)
  #  o$value <- as.matrix(o$value); colnames(o$value) <- colnames(model@y)


  return(list(par=o$par, value=values,fcount=fcount))
}

#include <cmath>

#include "misc.h"
#include "dt.h"

#include <RcppArmadillo.h>
#include <algorithm>
#include <math.h>
#include <limits.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


//' @title Rcpp implementation of Felzenszwalb distance transfom
//' @description Rcpp wrapper for the distance transform algorithm described in Felzenszwalb and Huttenlocher (2012)
//' @param x matrix of booleans of size \eqn{n x m} representing a (binary) image
//' @return A matrix of size \eqn{n x m} containing the distance transform result. Note that this function does not perform any checks on \code{x}.
//' @author Pedro Felzenszwalb for the header files \code{dt.h} and \code{misc.h} that do the work, Dario Azzimonti and Julien Bect for the wrapper.
//' @references  Felzenszwalb, P. F. and Huttenlocher, D. P. (2012). Distance Transforms of Sampled Functions. Theory of Computing, 8(19):415-428.
//' @examples
//' # Create an image with a square
//' nc = 256
//' nr = 256
//' xx = matrix(FALSE,ncol=nc,nrow=nr)
//' xx[(nr/16):(nr/16*15-1),nc/16]<-rep(TRUE,nr/16*14)
//' xx[(nr/16):(nr/16*15-1),nc/16*15]<-rep(TRUE,nr/16*14)
//' xx[nr/16,(nc/16):(nc/16*15-1)]<-rep(TRUE,nc/16*14)
//' xx[nr/16*15,(nc/16):(nc/16*15-1)]<-rep(TRUE,nc/16*14)
//' # Compute Distance transform
//' zz<- dtt_fast(xx)
//' \donttest{
//' # Plot the results
//' image(xx,col=grey.colors(20), main="Original image")
//' image(zz,col=grey.colors(20), main="Distance transform")
//' }
//' @export
// [[Rcpp::export]]
arma::mat dtt_fast(arma::mat x){
  int nrow = x.n_rows;
  int ncol = x.n_cols;

  // prepare the input for dt()
  image<uchar> *bw = new image<uchar>(ncol,nrow,false );

  for(int i=0; i<nrow; i++) // row index
    for(int j=0; j<ncol; j++) // column index
      imRef(bw,j,i) = (uchar)( x(i,j) );

  // compute the distance transform
  // (note: at this stage this is the *squared* distance)
  image<float> *y = dt(bw);

  // take the square root and stores the result into an arma object
  arma::mat z =  arma::zeros<arma::mat>(nrow,ncol);
  for(int i=0; i<nrow; i++) // row index
     for(int j=0; j<ncol; j++) // column index
       z(i,j) = sqrt( (double)imRef(y,j,i) );

  delete bw;
  delete y;
  return z;
}


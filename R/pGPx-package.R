#' @keywords internal
#' @details Generates posterior pseudo-realizations of Gaussian processes for excursion set estimation. The package provides posterior pseudo-realizations over large designs by simulating the field at few well chosen points and interpolating the result. The points are chosen minimizing the (posterior) expected distance in measure between the approximate excursion set and the full excursion set. The main functions in the package are: \describe{
#'    \item{\strong{Approximation:}}{ \itemize{
#'     \item \code{\link{optim_dist_measure}}: Given a \link[DiceKriging]{km} objects computes the optimal simulation points \eqn{e_1}, ... , \eqn{e_m} according to algorithm \code{A} or \code{B}.
#'     \item \code{\link{krig_weight_GPsimu}}: Given the simulations points and the interpolation points computes the kriging weights for the approximate process \eqn{\tilde{Z}} at the interpolation points.
#'     \item \code{\link{grad_kweights}}: Given the simulations points and the interpolation points returns the gradient of kriging weights with respect to the interpolation points.
#'     \item \code{\link{expDistMeasure}}: computes the expected distance in measure between the excursion set of the approximated process and the true excursion set.
#'    } }
#'    \item{\strong{Simulation:}}{ \itemize{
#'     \item \code{\link{simulate_and_interpolate}}: Generates \code{nsims} approximate posterior field realizations at \code{interpolatepoints} given the optimized simulation points.
#'    }  }
#'    \item{\strong{Applications:}}{ \itemize{
#'     \item \emph{Contour length}: the function \code{\link{compute_contourLength}} computes the excursion set contour length for each GP realization.
#'     \item \emph{Distance transform}: the function \code{\link{dtt_fast}} computes the distance transform of a binary image (Felzenszwalb and Huttenlocher, 2012) and the function \code{\link{DTV}} computes the distance transfom variability.
#'     \item \emph{Volumes}: the function \code{\link{computeVolumes}} computes the excursion volumes for each GP realization. It also applies a bias correction for approximate realizations.
#'    }  }
#' }
#' @note This work was supported in part by the Swiss National Science Foundation, grant numbers 146354, 167199 and the Hasler Foundation, grant number 16065. The author wishes to thank David Ginsbourger, Clément Chevalier and Julien Bect for the fruitful discussions.
#' @references Azzimonti, D., Bect, J., Chevalier, C., and Ginsbourger, D. (2016a). Quantifying uncertainties on excursion sets under a Gaussian random field prior. SIAM/ASA Journal on Uncertainty Quantification, 4(1):850–874.
#'
#' Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
#'
#' Azzimonti, D. and Ginsbourger, D. (2017). Estimating orthant probabilities of high dimensional Gaussian vectors with an application to set estimation. Journal of Computational and Graphical Statistics.
#'
#' Bolin, D. and Lindgren, F. (2015). Excursion and contour uncertainty regions for latent Gaussian models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 77(1):85--106.
#'
#' Chevalier, C. (2013). Fast uncertainty reduction strategies relying on Gaussian process models. PhD thesis, University of Bern.
#'
#'
#' Chevalier, C., Bect, J., Ginsbourger, D., Vazquez, E., Picheny, V., and Richet, Y. (2014). Fast kriging-based stepwise uncertainty reduction with application to the identification of an excursion set. Technometrics, 56(4):455–465.
#'
#' Felzenszwalb, P. F. and Huttenlocher, D. P. (2012). Distance Transforms of Sampled Functions. Theory of Computing, 8(19):415-428.
"_PACKAGE"

## usethis namespace: start
#' @import pbivnorm
#' @importFrom DiceKriging covMatrix covMat1Mat2 covVector.dx predict.km simulate
#' @importFrom KrigInv predict_nobias_km predict_update_km_parallel computeQuickKrigcov precomputeUpdateData vorob_threshold integration_design
#' @importFrom stats cov dist deriv pnorm rnorm runif model.matrix optim qnorm as.formula
#' @importFrom rgenoud genoud
#' @importFrom randtoolbox sobol
#' @importFrom pracma poly_length
#' @importFrom grDevices contourLines
#' @importFrom Rcpp evalCpp
#' @useDynLib pGPx
## usethis namespace: end
NULL

#' Full Kernel Estimation with Forward and Backward Lagged Covariates
#' 
#' A kernel weighting scheme to evaluate the effects of longitudinal covariates 
#'   on the occurrence of events when the time-dependent covariates are 
#'   measured intermittently. Regression parameter estimation uses full kernel 
#'   imputation of missing values with both forward and backward lagged 
#'   covariates.
#'
#' @param X An object of class data.frame. The structure of the data.frame must
#'  be \{patient ID, event time, event indicator\}. Patient IDs must be of class 
#'  integer or be able to be coerced to class integer without loss of 
#'  information. Missing values must be indicated as NA. The event indicator is
#'  1 if the event occurred; 0 if censored.
#' @param Z An object of class data.frame. The structure of the data.frame must 
#'   be \{patient ID, time of measurement, measurement(s)\}. Patient IDs must be 
#'   of class integer or be able to be coerced to class integer without loss of 
#'   information. Missing values must be indicated as NA.
#' @param tau An object of class numeric. The desired time point.
#' @param kType An object of class character indicating the type of smoothing 
#'   kernel to use in the estimating equation. Must be one of \{"epan", 
#'   "uniform", "gauss"\}, where "epan" is the Epanechnikov kernel and "gauss" 
#'   is the Gaussian kernel.
#' @param bw NULL or a numeric vector. If provided, the bandwidths for which 
#'   parameter estimates are to be obtained. If NULL, an optimal bandwidth will 
#'   be determined using an adaptive selection procedure. The range of the 
#'   bandwidth search space is taken to be 
#'   \eqn{2*(Q3 - Q1)*n^{-0.7}}{2*(Q3 - Q1)*n^{-0.7}} to 
#'   \eqn{2*(Q3 - Q1)*n^{-0.3}}{2*(Q3 - Q1)*n^{-0.3}}, where Q3 is the 0.75 
#'   quantile and Q1 is the 0.25 quantile of the measurement times for the 
#'   covariate and n is the effective number of patients, taken as the total 
#'   number of patients that experienced an event.
#' @param tol An object of class numeric. The minimum change in the regression 
#'   parameters deemed to indicate convergence of the Newton-Raphson method.
#' @param maxiter An object of class integer. The maximum number of iterations
#'  used to estimate regression parameters.
#' @param verbose An object of class logical. TRUE results in progress screen 
#'   prints.
#'
#' @returns A list is returned. If bandwidths are provided, each element is a 
#'   matrix, where the ith row corresponds to the ith bandwidth of input 
#'   argument \code{bw}, and the columns correspond to the model parameters. If 
#'   the bandwidth is determined internally, each element of the list is a 
#'   named vector calculated at the optimal bandwidth.
#'   \itemize{
#'     \item{betaHat }{The estimated model coefficients.}
#'     \item{stdErr  }{The standard error for each coefficient.}
#'     \item{zValue  }{The estimated z-value for each coefficient.}
#'     \item{pValue  }{The p-value for each coefficient.}
#'   }
#'
#'  If the bandwidth is determined internally, three additional list
#'  elements are returned:
#'   \itemize{
#'     \item{optBW }{The estimated optimal bandwidth.}
#'     \item{minMSE }{The mean squared error at the optimal bandwidth.}
#'     \item{MSE }{The vector of MSE for each bandwidth.}
#'   }
#'   
#' @references 
#'  Cao H., Churpek M. M., Zeng D., Fine J. P. 
#'  (2015).
#'  Analysis of the proportional hazards model with sparse longitudinal covariates.
#'  Journal of the American Statistical Association, 110, 1187-1196.
#' @seealso \code{\link{halfKernel}}, \code{\link{lastValue}}, \code{\link{nearValue}}
#'
#' @examples 
#'  data(SurvLongData)
#'
#'  exp <- fullKernel(X = X, Z = Z, tau = 1.0, bw = 0.015)
#'
#' @include kernelAuto.R kernelFixed.R
#' @export
fullKernel <- function(X, 
                       Z, 
                       tau,
                       kType = c("epan", "uniform", "gauss"), 
                       bw = NULL,
                       tol = 0.001,
                       maxiter = 100L, 
                       verbose = TRUE) {
  
  kType <- match.arg(kType)
  
  stopifnot(
    "`X` must be a data.frame with 3 columns" = !missing(X) && 
      is.data.frame(X) && ncol(X) == 3L,
    "`Z` must be a data.frame with at leat 3 columns" = !missing(Z) && 
      is.data.frame(Z) && ncol(Z) >= 3L,
    "`tau must be a scalar numeric" = !missing(tau) && is.numeric(tau) &&
      is.vector(tau) && length(tau) == 1L,
    "`bw` must be NULL or a numeric vector" = is.null(bw) ||
      {is.numeric(bw) && is.vector(bw)},
    "`tol` must be a positive scalar" = is.numeric(tol) && is.vector(tol) &&
      length(tol) == 1L && tol > 0.0,
    "`maxiter` must be an integer" = is.numeric(maxiter) && 
      isTRUE(all.equal(maxiter, round(maxiter))) && maxiter > 0,
    "`verbose` must be a logical" = is.logical(verbose)
  )

  if (is.null(bw)) {
    kernelAuto(X = X,
               Z = Z,
               tau = tau,
               kType = kType,
               tol = tol,
               maxiter = maxiter,
               scoreFunction = "scoreFull",
               verbose = verbose)
  } else {
    kernelFixed(X = X,
                Z = Z,
                tau = tau,
                bandwidth = bw,
                kType = kType,
                tol = tol,
                maxiter = maxiter,
                scoreFunction = "scoreFull",
                verbose = verbose)
  }
}

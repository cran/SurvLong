#' Half Kernel Estimation with Backward Lagged Covariates
#' 
#' A kernel weighting scheme to evaluate the effects of longitudinal covariates 
#'   on the occurrence of events when the time-dependent covariates are 
#'   measured intermittently. Regression parameter estimation using half kernel 
#'   imputation of missing values with backward lagged covariates. 
#'   
#' @inherit fullKernel params return references
#' @seealso \code{\link{fullKernel}}, \code{\link{lastValue}}, \code{\link{nearValue}}
#'
#' @examples 
#'  data(SurvLongData)
#'
#'  exp <- halfKernel(X = X, Z = Z, tau = 1.0, bw = 0.015)
#'
#' @include kernelAuto.R kernelFixed.R
#' @export
halfKernel <- function(X, 
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
               scoreFunction = "scoreHalf", 
               verbose = verbose)

  } else {

    kernelFixed(X = X,
                Z = Z,
                tau = tau,
                bandwidth = bw,
                kType = kType,
                tol = tol,
                maxiter = maxiter,
                scoreFunction = "scoreHalf", 
                verbose = verbose)

  }
}

#' Nearest Value Method
#' 
#' A simple approach to evaluate the effects of longitudinal covariates on the 
#'   occurrence of events when the time-dependent covariates are measured 
#'   intermittently. Regression parameters are estimated using the nearest 
#'   value to imputate missing values.
#'   
#' @inherit fullKernel params
#' @inherit lastValue return references
#' @seealso \code{\link{fullKernel}}, \code{\link{halfKernel}}, \code{\link{lastValue}}
#'
#' @examples
#'  data(SurvLongData)
#'  # A truncated dataset to keep example run time brief
#'  exp <- nearValue(X = X[1:100,], Z = Z, tau = 1.0)
#'
#' @importFrom stats pnorm
#' @include preprocessInputs.R betaEst.R scoreNVCF.R
#' @export
nearValue <- function(X, 
                      Z, 
                      tau,
                      tol = 0.001,
                      maxiter = 100L, 
                      verbose = TRUE) {

  # Process and verify input datasets                                        #
  pre <- preprocessInputs(data.x = X, data.z = Z)

  X <- pre$data.x
  Z <- pre$data.z
  nCov <- ncol(Z) - 2L

  rm(pre)

  # Calculate parameter estimates and standard deviations                    #
  bHat <- betaEst(Z = Z,
                  X = X, 
                  tau = tau,
                  tol = tol,
                  h = 0,
                  kType = "epan",
                  betaGuess = NULL,
                  maxiter = maxiter,
                  scoreFunction = "scoreNVCF")

  score <- scoreNVCF(beta = bHat,
                     Z = Z,
                     X = X, 
                     tau = tau,
                     h = 0,
                     kType = "epan")

  invdU <- tryCatch(solve(score$dUdBeta), 
                    error = function(e) {
                      stop("unable to invert derivative of estimating equation\n\t",
                           e$message, call. = FALSE)
                    })

  sig <- invdU %*% score$mMatrix %*% invdU

  sdVec <- sqrt(diag(sig)) 

  if (verbose) {
    # Generate results matrix
    results <- matrix(0.0, nrow = nCov, ncol = 4L,
                      dimnames = list(paste0("beta",1L:{nCov}),
                                      c("estimate","stdErr","z-value","p-value")))
    results[, 1L] <- bHat
    results[, 2L] <- sdVec
    results[, 3L] <- bHat / sdVec
    results[, 4L] <- 2.0 * stats::pnorm(-abs(results[, 3L]))
    
    print(results)
  }

  zv <- bHat / sdVec
  pv <- 2.0 * stats::pnorm(-abs(zv))

  list("betaHat" = matrix(bHat, nrow = 1L),
       "stdErr"  = sdVec,
       "zValue" = zv,
       "pValue" = pv)

}

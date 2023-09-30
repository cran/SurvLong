#' Last Value Carried Forward Method
#'
#' A simple approach to evaluate the effects of longitudinal covariates on the 
#'   occurrence of events when the time-dependent covariates are measured 
#'   intermittently. Regression parameter are estimated using last value 
#'   carried forward imputation of missing values.
#'   
#' @inherit fullKernel params references
#' @returns A list 
#'   \itemize{
#'     \item{betaHat }{The estimated model coefficients.}
#'     \item{stdErr  }{The standard error for each coefficient.}
#'     \item{zValue  }{The estimated z-value for each coefficient.}
#'     \item{pValue  }{The p-value for each coefficient.}
#'   }
#'   
#'  
#' @seealso \code{\link{fullKernel}}, \code{\link{halfKernel}}, \code{\link{nearValue}}
#' 
#' @examples 
#'  data(SurvLongData)
#'  # A truncated dataset to keep example run time brief
#'  exp <- lastValue(X = X[1:200,], Z = Z, tau = 1.0)
#'  
#' @importFrom stats pnorm
#' @include preprocessInputs.R betaEst.R scoreLVCF.R
#' @export
lastValue <- function(X, 
                      Z, 
                      tau,
                      tol = 0.001,
                      maxiter = 100L, 
                      verbose = TRUE) {

  # Process and verify input datasets
  pre <- preprocessInputs(data.x = X, data.z = Z)

  X <- pre$data.x
  Z <- pre$data.z
  nCov <- ncol(Z) - 2L

  rm(pre)

  # Calculate parameter estimates and standard deviations
  bHat <- betaEst(Z = Z,
                  X = X, 
                  tau = tau,
                  h = 0,
                  kType = "epan",
                  tol = tol,
                  betaGuess = NULL,
                  maxiter = maxiter,
                  scoreFunction = "scoreLVCF")

  score <- scoreLVCF(beta = bHat,
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
       "stdErr" = sdVec,
       "zValue" = zv,
       "pValue" = pv )

}

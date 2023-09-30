#' Fixed bandwidth algorithm for full and half kernel methods.
#'
#' @noRd
#' @param Z An object of class data.frame. The structure of the data.frame must 
#'   be \{patient ID, time of measurement, measurement(s)\}. Patient IDs must be 
#'   of class integer or be able to be coerced to class integer without loss of 
#'   information. Missing values must be indicated as NA.
#' @param X An object of class data.frame. The structure of the data.frame must 
#'   be \{patient ID, event time, event indicator\}. Patient IDs must be of 
#'   class integer or be able to be coerced to class integer without loss of 
#'   information. Missing values must be indicated as NA.
#' @param tau An object of class numeric. The desired time point.
#' @param bandwidth An object of class numeric. The bandwidth value(s) at which 
#'   parameters are to be estimated.
#' @param kType An object of class character indicating the type of smoothing 
#'   kernel to use in the estimating equation. Must be one of \{"epan", 
#'   "uniform", "gauss"\}, where "epan" is the Epanechnikov kernel and "gauss" 
#'   is the Gaussian kernel.
#' @param tol An object of class numeric. The maximum allowed change in 
#'   parameter estimates, beyond which the parameter estimates are deemed to 
#'   have converged.
#' @param maxiter An object of class numeric. The maximum number of iterations 
#'   allowed to attain convergence.
#' @param scoreFunction An object of class character. The name of the function 
#'   to be used to calculate the score.
#'   
#' @returns A list
#'   \itemize{
#'     \item{betaHat }{The estimated model coefficients.}
#'     \item{stdErr }{The standard error for each coefficient.}
#'     \item{zValue }{The estimated z-value for each coefficient.}
#'     \item{pValue }{The p-value for each coefficient.}
#'  }
#'  
#' @include betaEst.R preprocessInputs.R 
#' @include scoreFull.R scoreHalf.R scoreLVCF.R scoreNVCF.R
#' @importFrom stats pnorm
#' @keywords internal
kernelFixed <- function(X, 
                        Z, 
                        tau, 
                        bandwidth,
                        kType,
                        tol,
                        maxiter,
                        scoreFunction, 
                        verbose){

  # Process and verify input datasets
  pre <- preprocessInputs(data.x = X, data.z = Z)

  X <- pre$data.x
  Z <- pre$data.z
  nCov <- ncol(Z) - 2L

  rm(pre)

  # Determine number of bandwidths provided by user
  lbd <- length(bandwidth)

  # initialize matrices for parameter estimates and standard deviations
  bHat <- matrix(0.0, nrow = lbd, ncol = nCov,
                 dimnames = list(NULL,colnames(Z)[3L:(nCov+2L)]))

  sdVec <- bHat

  guess <- NULL

  results <- matrix(0.0, nrow = nCov, ncol = 4L,
                    dimnames = list(paste0("beta", 0L:{nCov-1L}),
                                    c("estimate", "stdErr", "z-value", "p-value")))
  for (bd in 1L:lbd) {

    if (verbose) message("Bandwidth: ", bandwidth[bd])

    bHat[bd, ] <- betaEst(Z = Z,
                          X = X, 
                          tau = tau, 
                          h = bandwidth[bd],
                          kType = kType,
                          betaGuess = guess,
                          tol = tol,
                          maxiter = maxiter,
                          scoreFunction = scoreFunction)

    guess <- bHat[bd, ]

    argList <- list("beta" = bHat[bd, ],
                    "Z" = Z,
                    "X" = X, 
                    "tau" = tau, 
                    "h" = bandwidth[bd],
                    "kType" = kType)

    score <- do.call(scoreFunction, args = argList)

    invdU <- tryCatch(solve(score$dUdBeta), 
                      error = function(e) {
                        stop("Unable to invert derivative of estimating equation.\n\t",
                             e$message, call. = FALSE)
                      })

    sig <- invdU %*% score$mMatrix %*% invdU

    sdVec[bd, ] <- sqrt(diag(sig)) 

    if (verbose) {
      results[ ,1L] <- bHat[bd, ]
      results[ ,2L] <- sdVec[bd, ]
      results[ ,3L] <- bHat[bd, ] / sdVec[bd, ]
      results[ ,4L] <- 2.0 * stats::pnorm(-abs(results[, 3L]))
      
      print(results)
      cat("\n")
    }

  }

  list("betaHat" = bHat,
       "stdErr"  = sdVec,
       "zValue" =  bHat / sdVec,
       "pValue" = 2.0 * stats::pnorm(-abs(bHat / sdVec)))

}

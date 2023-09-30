#' Auto Tune algorithm for full and half kernel methods.
#'
#' @noRd
#' @param Z An object of class data.frame. The structure of the data.frame must 
#'   be \{patient ID, time of measurement, measurement(s)\}. Patient IDs must 
#'   be of class integer or be able to be coerced to class integer without loss 
#'   of information. Missing values must be indicated as NA.
#' @param X An object of class data.frame. The structure of the data.frame must 
#'   be \{patient ID, event time, event indicator\}. Patient IDs must be of class 
#'   integer or be able to be coerced to class integer without loss of 
#'   information. Missing values must be indicated as NA.
#' @param tau An object of class numeric. The desired time point.
#' @param kType An object of class character indicating the type of smoothing 
#'   kernel to use in the estimating equation. Must be one of \{"epan", 
#'   "uniform", "gauss"\}, where "epan" is the Epanechnikov kernel and "gauss" 
#'   is the Gaussian kernel.
#' @param tol An object of class numeric. The maximum allowed change in 
#' parameter estimates, beyond which the parameter estimates are deemed to have 
#' converged.
#' @param maxiter An object of class numeric. The maximum number of iterations 
#' allowed to attain convergence.
#' @param scoreFunction An object of class character. The name of the function
#'  to be used to calculate the score.
#' @returns a list
#'   \itemize{
#'     \item{betaHat }{The estimated model coefficients.}
#'     \item{stdErr }{The standard error for each coefficient.}
#'     \item{zValue }{The estimated z-value for each coefficient.}
#'     \item{pValue }{The p-value for each coefficient.}
#'     }
#' If the bandwidth is determined automatically, two additional list 
#' elements are returned:
#'   \itemize{
#'     \item{optBW }{The estimated optimal bandwidth.}
#'     \item{minMSE }{The mean squared error at the optimal bandwidth.}
#'     }
#'     
#' @include betaEst.R preprocessInputs.R
#' @include scoreFull.R scoreHalf.R scoreLVCF.R scoreNVCF.R
#' @importFrom stats lm pnorm quantile
#' @keywords internal
kernelAuto <- function(X, 
                       Z, 
                       tau, 
                       kType,
                       tol,
                       maxiter,
                       scoreFunction, 
                       verbose) {

  # Process and verify input datasets
  pre <- preprocessInputs(data.x = X, data.z = Z)

  X <- pre$data.x
  Z <- pre$data.z
  nCov <- ncol(Z) - 2L

  rm(pre)

  # Determine search space for bandwidth
  effn <- sum(X[, 3L])
  limit <- 2.0 * (stats::quantile(Z[, 2L], 0.75) - stats::quantile(Z[, 2L], 0.25))
  bw <- seq(from = limit * (effn)^(-0.7), 
            to = limit * (effn)^(-0.3),
            length.out = 50)
  lbd <- length(bw)

  # initialize matrix for parameter estimates
  betaHat0 <- matrix(0.0, nrow = lbd, ncol = nCov)

  # initialize matrices for variance estimate
  hatV <- matrix(0.0, nrow = lbd, ncol = nCov)

  # Randomly assign each patient to group 1 or 2
  cvlabel <- sample(1L:2L, nrow(X), replace = TRUE)

  # Estimate parameters at each bandwidth

  guess0 <- NULL
  guess1 <- NULL
  guess2 <- NULL

  for( bd in 1L:lbd ) {

    betaHat0[bd, ] <- betaEst(Z = Z,
                              X = X, 
                              tau = tau, 
                              h = bw[bd],
                              kType = kType,
                              betaGuess = guess0,
                              tol = tol,
                              maxiter = maxiter,
                              scoreFunction = scoreFunction)

    guess0 <- betaHat0[bd, ]

    tst <- cvlabel == 1L

    betaHat1 <- betaEst(Z = Z,
                        X = X[tst, , drop = FALSE], 
                        tau = tau, 
                        h = bw[bd],
                        kType = kType,
                        betaGuess = guess1,
                        tol = tol,
                        maxiter = maxiter,
                        scoreFunction = scoreFunction)

    guess1 <- betaHat1

    betaHat2 <- betaEst(Z = Z,
                        X = X[!tst, , drop = FALSE], 
                        tau = tau, 
                        h = bw[bd],
                        kType = kType,
                        betaGuess = guess2,
                        tol = tol,
                        maxiter = maxiter,
                        scoreFunction = scoreFunction)

    guess2 <- betaHat2

    betaDiff <- betaHat2 - betaHat1

    hatV[bd, ] <- nrow(X) * bw[bd] * (betaDiff * betaDiff) * 0.25

  }

  # Estimate slope of bias expression; Calculate MSE
  hatC <- numeric(nCov)
  MSE <- numeric(lbd)

  if (scoreFunction == "scoreHalf") {

    for (p in 1L:nCov) {
      hatC[p] <- stats::lm( betaHat0[, p] ~ bw )$coef[2L]
    }

    # MSE = hat(C) * hat(C) * bw^2 + hat(V)
    MSE <- bw^2 * drop(hatC %*% hatC) + rowSums(hatV)

  } else if (scoreFunction == "scoreFull") {
    # 9/30/2023 This was previously implemented as y ~ bw^2,  which is not correct
    bw2 <- bw^2
    for (p in 1L:nCov) {
      hatC[p] <- stats::lm( betaHat0[, p] ~ bw2 )$coef[2L]
    }

    # MSE = hat(C) * hat(C) * bw^4 + hat(V)
    MSE <- bw^4 * drop(hatC %*% hatC) + rowSums(hatV)

  }

  # Identify minimum MSE
  tst <- MSE > 0.0
  if (!any(tst)) stop("no positive MSE values", call. = FALSE)
  MSE[!tst] <- NA_real_

  opt_h <- which.min(MSE)

  if (length(opt_h) > 1L) {
    warning("Multiple minimums. Smallest bandwidth used.", call. = FALSE)
    opt_h <- opt_h[1L]
  }

  if (isTRUE(all.equal(opt_h, bw[1L])) || isTRUE(all.equal(opt_h, bw[lbd]))) {
    warning("Minimum is at bandwidth boundary.", call. = FALSE)
  }

  minMSE <-  MSE[opt_h]

  # Estimate parameters and standard error at optimal bandwidth
  bHat <- betaEst(Z = Z,
                  X = X, 
                  tau = tau, 
                  h = bw[opt_h],
                  kType = kType,
                  betaGuess = guess0,
                  tol = tol,
                  maxiter = maxiter,
                  scoreFunction = scoreFunction)

  names(bHat) <- colnames(Z)[3L:(2L + nCov)]

  argList <- list("beta" = bHat,
                  "Z" = Z,
                  "X" = X, 
                  "tau" = tau, 
                  "h" = bw[opt_h],
                  "kType" = kType)

  score <- do.call(scoreFunction, args = argList)

  invdU <- tryCatch(solve(score$dUdBeta), 
                    error = function(e) {
                      stop("Unable to invert derivative of estimating equation.\n\t",
                           e$message, call. = FALSE)
                    })

  sig <- invdU %*% score$mMatrix %*% invdU

  sdVec <- sqrt(diag(sig)) 

  names(sdVec) <- names(bHat)

  # Generate results matrix
  
  if (verbose) {
    results <- cbind("estimate" = bHat,
                     "stdErr" = sdVec,
                     "z-value" = bHat / sdVec,
                     "p-value" = 2.0 * stats::pnorm(-abs(bHat / sdVec)))
    rownames(results) <- paste0("beta", 0L:{nCov - 1L})
    
    message("Bandwidth search range: ", bw[1L], " - ", bw[lbd])
    message("Optimal bandwidth: ", bw[opt_h])
    message("MSE: ", minMSE)
    print(results)
  }

  zv <- bHat / sdVec
  pv <- 2.0 * stats::pnorm(-abs(zv))
  MSE <- matrix(MSE, ncol = 1L)
  rownames(MSE) <- format(bw, digits = 4L)

  list("betaHat" = bHat,
       "stdErr" = sdVec,
       "zValue" = zv,
       "pValue" = pv,
       "optBW" = bw[opt_h],
       "minMSE" =  minMSE,
       "MSE" = MSE)

}

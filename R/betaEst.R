#' Newton Raphson algorithm for all methods.
#' 
#' @noRd
#' @param Z An object of class data.frame. The structure of the data.frame must                    #
#'  be \{patient ID, time of measurement, measurement(s)\}. Patient IDs must be 
#'  of class integer or be able to be coerced to class integer without loss of 
#'  information. Missing values must be indicated as NA.
#' @param X An object of class data.frame. The structure of the data.frame must 
#'   be \{patient ID, event time, event indicator\}. Patient IDs must be of 
#'   class integer or be able to be coerced to class integer without loss of 
#'   information. Missing values must be indicated as NA.
#' @param tau An object of class numeric. The desired time point.
#' @param h An object of class numeric. The bandwidth.
#' @param kType An object of class character indicating the type of smoothing 
#'   kernel to use in the estimating equation. Must be one of \{"epan", 
#'   "uniform", "gauss"\}, where  "epan" is the Epanechnikov kernel and "gauss" 
#'   is the Gaussian kernel.
#' @param betaGuess An object of class numeric or NULL. If numeric, beta will 
#'   be initialized to the values provided.
#' @param tol An object of class numeric. The maximum allowed change in 
#'   parameter estimates, beyond which the parameter estimates are deemed to 
#'   have converged.
#' @param maxiter An object of class numeric. The maximum number of iterations 
#'   allowed to attain convergence.
#' @param scoreFunction An object of class character. The name of the function 
#'   to be used to calculate the score.
#'  
#' @returns An numeric vector containing the parameter estimates.
#' 
#' @keywords internal
betaEst <- function(Z,  
                    X,  
                    tau,  
                    h,
                    kType,
                    betaGuess,
                    tol,
                    maxiter,
                    scoreFunction) {


  # If a starting value for Newton-Raphson provided, use. Else, initialize
  # to small positive value (0.01).
  if (is.null(betaGuess)) {
    beta <- rep(0.01, ncol(Z) - 2L)
  } else {
    beta <- betaGuess
  }

  iter <- 0L

  while (TRUE) {

    # Calculate Score Function
    argList <- list("beta" = beta, 
                    "Z" = Z, 
                    "X" = X, 
                    "tau" = tau, 
                    "h" = h, 
                    "kType" = kType)
    
    Lvec <- do.call(scoreFunction, argList)
    
    change <- tryCatch(solve(Lvec$dUdBeta, Lvec$U),
                        error = function(e) {
                          stop("Unable to invert matrix\n\t",
                               e$message, call. = FALSE)
                        })
    
    if (any(is.na(change))) {
      stop("NAs encountered in Newton-Raphson", call. = FALSE)
    }

    beta.hat <- beta - change

    # Determine if parameter estimates have converged.
    test <- abs(change / beta)
    test[abs(beta) < 0.001] <- abs(change)[abs(beta) < 0.001]
    if (all(test < tol)) break

    beta <- beta.hat

    # Increment iterations and verify that maximum not yet reached
    iter <- iter + 1L
    if (iter >= maxiter) {
      warning(paste("Parameter estimates did not converge within", 
                    maxiter, "iterations."), call. = FALSE)
      break
    }
  }

  beta
}

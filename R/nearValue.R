#******************************************************************************#
# Public calling function for nearest value method                             #
#******************************************************************************#
#                                                                              #
# Inputs                                                                       #
#                                                                              #
#  X              an object of class data.frame.                               #
#                 The structure of the data.frame must be                      #
#                 \{patient ID, event time, event indicator\}.                 #
#                 Patient IDs must be of class integer or be able to be        #
#                 coerced to class integer without loss of information.        #
#                 Missing values must be indicated as NA.                      #
#                                                                              #
#  Z              an object of class data.frame.                               #
#                 The structure of the data.frame must be                      #
#                 \{patient ID, time of measurement, measurement(s)\}.         #
#                 Patient IDs must be of class integer or be able to be        #
#                 coerced to class integer without loss of information.        #
#                 Missing values must be indicated as NA.                      #
#                                                                              #
#  tau            an object of class numeric.                                  #
#                 The desired time point.                                      #
#                                                                              #
#  tol            an object of class numeric.                                  #
#                 maximum allowed change in parameter estimates, beyond which  #
#                 the parameter estimates are deemed to have converged.        #
#                                                                              #
#  maxiter        an object of class numeric.                                  #
#                 maximum number of iterations allowed to attain convergence   #
#                                                                              #
#  verbose        an object of class logical.                                  #
#                 TRUE indicates that screen prints are generated              #
#                                                                              #
#  Outputs                                                                     #
#                                                                              #
#  Returns a list                                                              #
#                                                                              #
# betaHat The estimated model coefficients.                                    #
# stdErr  The standard error for each coefficient.                             #
# zValue  The estimated z-value for each coefficient.                          #
# pValue  The p-value for each coefficient.                                    #
#                                                                              #
#******************************************************************************#
nearValue <- function(X, 
                      Z, 
                      tau,
                      tol = 0.001,
                      maxiter = 100, 
                      verbose = TRUE){

  #--------------------------------------------------------------------------#
  # Process and verify input datasets                                        #
  #--------------------------------------------------------------------------#
  pre <- preprocessInputs(data.x = X, data.z = Z)

  X <- pre$data.x
  Z <- pre$data.z
  nCov <- ncol(Z) - 2L

  rm(pre)

  #--------------------------------------------------------------------------#
  # Calculate parameter estimates and standard deviations                    #
  #--------------------------------------------------------------------------#
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

  invdU <- try(solve(score$dUdBeta), silent = TRUE)

  if( is(invdU, 'try-error') ) {
    cat("unable to invert derivative of estimating equation\n")
    stop(attr(invdU,"condition"))
  }

  sig <- invdU %*% (score$mMatrix) %*% invdU

  sdVec <- sqrt(diag(sig)) 

  #--------------------------------------------------------------------------#
  # Generate results matrix                                                  #
  #--------------------------------------------------------------------------#
  results <- matrix(data = 0.0,
                    nrow = nCov,
                    ncol = 4L,
                    dimnames = list(paste("beta",1L:{nCov},sep=""),
                                    c("estimate","stdErr","z-value","p-value")))

  results[,1L] <- bHat
  results[,2L] <- sdVec
  results[,3L] <- bHat/sdVec
  results[,4L] <- 2.0*pnorm(-abs(results[,3L]))

  if (verbose) print(results)

  zv <- bHat/sdVec
  pv <- 2.0*pnorm(-abs(zv))

  return( list( "betaHat" = matrix(bHat,nrow=1L),
                "stdErr"  = sdVec,
                "zValue" = bHat/sdVec,
                "pValue" = pv ) )

}

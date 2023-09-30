#' Score calculation for half kernel method
#'
#' @noRd
#' @param beta An object of class numeric. The parameter estimate(s).
#' @param Z An object of class data.frame. The structure of the data.frame must 
#'   be \{patient ID, time of measurement, measurement(s)\}. Missing values 
#'   should have been previously set to 0 or removed.
#' @param X An object of class data.frame. The structure of the data.frame must
#'   be \{patient ID, event time, event indicator\}. Missing values 
#'   should have been previously set to 0 or removed.
#' @param tau An object of class numeric. The desired time point.
#' @param h An object of class numeric. The  kernel bandwidth.
#' @param kType An object of class character indicating the type of smoothing 
#'   kernel to use in the estimating equation. Must be one of \{"epan", 
#'   "uniform", "gauss"\}, where "epan" is the Epanechnikov kernel and "gauss" 
#'   is the Gaussian kernel.
#' @param ... Ignored.
#' 
#' @returns A list 
#'  \itemize {
#'    \item{U }{An object of class numeric. The Score function(s).}
#'    \item{dUdBeta }{An object of class numeric matrix. The derivative of the 
#'      Score function.}
#'    \item{mMatrix }{An object of class numeric matrix. Sigma}
#'  }
#'    
#' @include local_kernel.R
#' @keywords internal
scoreHalf <- function(beta, 
                      Z, 
                      X,  
                      tau,  
                      h,
                      kType, ...) {

  p <- ncol(Z) - 2L

  Zp <- data.matrix(Z[, -c(1L:2L), drop = FALSE])
  Zid <- Z[, 1L]
  Ztime <- Z[, 2L]
  rm(Z)

  dUdBeta <- matrix(0.0, nrow = p, ncol = p)
  Mmatrix <- matrix(0.0, nrow = p, ncol = p)
  Uvec <- numeric(p)

  n <- nrow(X)

  # X^T X for each participant
  outerOnce <- apply(Zp, 1L, tcrossprod) |> t()
  if (ncol(outerOnce) == nrow(Zp) && nrow(outerOnce) == p * p) {
    outerOnce <- t(outerOnce)
  }
  
  for (i in 1L:n) {

    # If the time is censored, do not include in summation
    if( X[i, 3L] < 0.5 ) next

    # If the time is greater than the integration limit, do not include in
    # the summation
    time <- X[i, 2L]

    if (time > tau) next

    # Keep only those covariates with measurement time less equal time
    use <- (Ztime <= time) & Zid == X[i, 1L]

    if (!any(use)) next

    # Calculate the S Function
    kern <- local_kernel(t = {time - Ztime}, 
                         h = h,  
                         kType = kType)

    # Identify patients still at risk (t >= time)
    ptIDs <- X[time <= X[, 2L], 1L]

    # Identify the covariates for this subset of patients
    ZptIDs <- Zid %in% ptIDs & (Ztime <= time)

    cova <- Zp[ZptIDs, , drop = FALSE]

    prod <- kern[ZptIDs] * exp(cova %*% beta)

    s0 <- sum(prod)

    if( (s0 > -1.5e-8) && (s0 < 1.5e-8) ) next

    s1 <- colSums(prod[, 1L] * cova)

    s2 <- matrix(data = colSums(outerOnce[ZptIDs, , drop = FALSE] * prod[, 1L]),
                 nrow = p,
                 ncol = p)

    # Calculate U and dUdBeta                                              #
    ZmZ <- sweep(x = Zp[use, , drop = FALSE], 
                 MARGIN = 2L,
                 STATS = s1 / s0,
                 FUN = "-")

    tmp <- colSums(kern[use] * ZmZ)

    Mmatrix <- Mmatrix + tcrossprod(tmp)

    Uvec <- Uvec + tmp

    dUdBeta <- dUdBeta + sum(kern[use]) * (tcrossprod(s1) - s2 * s0) / (s0 * s0)

  }

  list("U" = Uvec, "dUdBeta" = dUdBeta, "mMatrix" = Mmatrix)

}


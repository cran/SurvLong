#' Score calculation for last value carried forward method
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
scoreLVCF <- function(beta, 
                      Z,  
                      X,  
                      tau, ...) {

  p <- ncol(Z) - 2L

  Lmat <- matrix(0.0, nrow = p, ncol = p)
  Mmatrix <- matrix(0.0, nrow = p, ncol = p)
  Uvec <- numeric(p)

  n <- nrow(X)
  nZ <- nrow(Z)
  
  extract_func <- function(x) {
    it <- which.max(x[, 2L])
    x[it, ]
    }

  for (i in 1L:n) {

    # If the time is censored, do not include in summation
    if (X[i, 3L] < 0.5) next

    # If the time is greater than the integration limit, do not include in
    # the summation
    time <- X[i, 2L]

    if (time > tau) next

    # Calculate the S Function

    # Identify patients still at risk (t >= time)
    ptIDs <- X[time <= X[, 2L], 1L]

    # Identify the covariates for this subset of patients
    ZptIDs <- Z[, 1L] %in% ptIDs & (Z[, 2L] <= time)

    if (!any(ZptIDs)) next

    Z2 <- Z[ZptIDs, , drop = FALSE]

    # Identify largest times and accept those covariates
    cova <- by(data = Z2,
               INDICES = Z2[, 1L],
               FUN = extract_func,
               simplify = FALSE) |> 
      unlist() |> 
      matrix(ncol = {p + 2L}, byrow = TRUE)

    IDs <- cova[, 1L] == X[i, 1L]
    if (sum(IDs) != 1L) next
    cova <- cova[, -c(1L:2L), drop = FALSE]

    prod <- exp(cova %*% beta)

    s0 <- sum(prod)

    if ((s0 > -1.5e-8) && (s0 < 1.5e-8)) next

    s1 <- colSums(prod[, 1L] * cova)

    tst <- matrix(apply(cova, 1L, tcrossprod), nrow = p * p)
    Zp <- rowSums(sweep(x = tst,
                        MARGIN = 2L,
                        STATS = prod[, 1L],
                        FUN = "*"))
    s2 <- matrix(Zp, nrow = p, ncol = p)

    tmp <- cova[IDs,] - s1 / s0

    Mmatrix <- Mmatrix + tcrossprod(tmp)

    Uvec <- Uvec + tmp

    Lmat <- Lmat + (tcrossprod(s1) - s2 * s0) / (s0 * s0)

  }

  list("U" = Uvec, "dUdBeta" = Lmat, "mMatrix" = Mmatrix)
}
#' Verify and pre-process inputs
#' 
#' @noRd
#' @param data.x An object of class data.frame. The structure of the data.frame 
#'   must be \{patient ID, event time, event indicator\}. Missing values must 
#'   be indicated as NA.
#' @param data.z An object of class data.frame. The structure of the data.frame 
#'   must be \{patient ID, time of measurement, measurement(s)\}. Missing 
#'   values must be indicated as NA.
#' @returns A list containing 
#'   \itemize{
#'     \item{data.x }{Same as input with: NAs responses removed}
#'     \item{data.z }{Same as input with: missing data cases set to 0}
#'     }
#'   
#' @keywords internal
preprocessInputs <- function(data.x, data.z) {

  stopifnot(
    "`X` must include {ID, time, delta}" = is.data.frame(data.x) &&
      ncol(data.x) == 3L,
    "`Z` must include {ID, time, measurement(s)}" = is.data.frame(data.z) &&
      ncol(data.z) >= 3L
  )
  
  if( !is.integer(data.z[, 1L]) ) {
    data.z[,1L] <- as.integer(round(data.z[, 1L], 0))
    message("Patient IDs in `Z` were coerced to integer.\n")
  }
  
  if( !is.integer(data.x[,1L]) ) {
    data.x[, 1L] <- as.integer(round(data.x[, 1L], 0))
    message("Patient IDs in `X` were coerced to integer.\n")
  }

  # Remove any cases for which all covariates are NA
  rmRow <- apply(data.z, 1L, function(x) { all(is.na(x)) })
  data.z <- data.z[!rmRow, ]

  # Set missing cases to 0.0
  data.z[is.na(data.z)] <- 0.0

  if (any(data.z[, 2L] < {-1.5e-8})) {
    stop("time variable is negative in `Z`", call. = FALSE)
  }
  
  # Remove any cases for which response is NA
  tst <- is.na(data.x[, 2L])
  data.x <- data.x[!tst, ]

  if (any(data.x[, 2L] < {-1.5e-8})) {
    stop("time is negative in `X`", call. = FALSE)
  }

  list("data.x" = data.x, "data.z" = data.z)
}

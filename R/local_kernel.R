# Kernel functions

#' @noRd
#' @param t An object of class numeric. The time point(s) at which the kernel is 
#'   to be calculated.
#' @param  h An object of class  numeric. The kernel bandwidth.
#' @param kType An object of class character indicating the type of smoothing 
#'   kernel to use in the estimating equation. Must be one of \{"epan", 
#'   "uniform", "gauss"\}, where "epan" is the Epanechnikov kernel and "gauss" 
#'   is the Gaussian kernel.
#'   
#' @returns An object of class numeric.
#' 
#' @keywords internal
local_kernel <- function(t, h, kType) {
  
  switch(kType,
         "epan" = .epanechnikov(t / h) / h,
         "uniform" = .uniform(t / h) / h,
         "gauss" = .gauss(t / h) / h,
         stop("unsupported kernel", call. = FALSE))
}

#' Epanechnikov Kernel
#' @noRd
#' @keywords internal
.epanechnikov <- function(t) {

  tst <- (-1.0 <= t) & (t <= 1.0 )

  kt <- 0.75 * (1.0 - t * t)
  kt[!tst] <- 0.0

  kt
}

#' Uniform Kernel
#' @noRd
#' @keywords internal
.uniform <- function(t) {

  tst <- (-1.0 <= t) & (t <= 1.0 )
  kt <- t
  kt[tst] <- 0.5
  kt[!tst] <- 0.0

  kt
}

#' Gaussian Kernel
#' @noRd
#' @keywords internal
.gauss <- function(t) {
  exp(-t * t * 0.5) / sqrt(2.0 * pi)
}

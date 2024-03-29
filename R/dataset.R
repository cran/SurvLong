#' Generated Sparse Longitudinal Data
#' 
#' For the purposes of the package examples, the dataset was adapted from the 
#'   numerical simulations of the original manuscript. 
#' 
#' Data was generated for 400 subjects. The total number of covariate observation 
#'   times was Poisson distributed with intensity rate 8. The covariate 
#'   observation times are generated from a uniform distribution Unif(0,1) 
#'   independently. The covariate process is piecewise constant, with values 
#'   being multivariate normal with mean 0, variance 1 and correlation 
#'   \eqn{\exp(-|i - j|/20)}{exp(-|i - j|/20)}. The survival time were generated 
#'   from the Cox model 
#'   \eqn{\lambda(t | Z(r), r \le t) = \lambda_0 \exp(\beta Z(t))}{
#'   lambda{t|Z(r),r<=t}=lambda0 exp(beta Z(t))}, where \eqn{\beta}{beta} = 1.5, 
#'   and \eqn{\lambda_0}{lambda0} = 1.0. Covariates are dataset Z. Event times 
#'   and indicators are dataset X.
#'   
#' @name SurvLongData
#' @aliases X Z
#'   
#' @format 
#'  X is a data frame with 400 observations on the following 3 variables.
#'  \describe{
#'    \item{\code{ID}}{patient identifier, there are 400 patients.}
#'    \item{\code{Time}}{the time to event or censoring}
#'    \item{\code{Delta}}{a numeric vector with 0 denoting censoring and 1 event}
#'  }
#'  Z is a data frame with 3237 observations on the following 3 variables.
#'  \describe{
#'    \item{\code{ID}}{patient identifier, there are 400 patients.}
#'    \item{\code{obsTime}}{the covariate observation times.}
#'    \item{\code{x1}}{the covariate generated through a piecewise constant function.}
#'  }
#'  
#' @references 
#'  Cao H., Churpek M. M., Zeng D., Fine J. P. 
#'  (2015).
#'  Analysis of the proportional hazards model with sparse longitudinal covariates.
#'  Journal of the American Statistical Association, 110, 1187-1196.
# @keywords dataset
NULL
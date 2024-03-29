% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nearValue.R
\name{nearValue}
\alias{nearValue}
\title{Nearest Value Method}
\usage{
nearValue(X, Z, tau, tol = 0.001, maxiter = 100L, verbose = TRUE)
}
\arguments{
\item{X}{An object of class data.frame. The structure of the data.frame must
be \{patient ID, event time, event indicator\}. Patient IDs must be of class 
integer or be able to be coerced to class integer without loss of 
information. Missing values must be indicated as NA. The event indicator is
1 if the event occurred; 0 if censored.}

\item{Z}{An object of class data.frame. The structure of the data.frame must 
be \{patient ID, time of measurement, measurement(s)\}. Patient IDs must be 
of class integer or be able to be coerced to class integer without loss of 
information. Missing values must be indicated as NA.}

\item{tau}{An object of class numeric. The desired time point.}

\item{tol}{An object of class numeric. The minimum change in the regression 
parameters deemed to indicate convergence of the Newton-Raphson method.}

\item{maxiter}{An object of class integer. The maximum number of iterations
used to estimate regression parameters.}

\item{verbose}{An object of class logical. TRUE results in progress screen 
prints.}
}
\value{
A list 
  \itemize{
    \item betaHat: The estimated model coefficients.
    \item stdErr: The standard error for each coefficient.
    \item zValue: The estimated z-value for each coefficient.
    \item pValue: The p-value for each coefficient.
  }
}
\description{
A simple approach to evaluate the effects of longitudinal covariates on the 
  occurrence of events when the time-dependent covariates are measured 
  intermittently. Regression parameters are estimated using the nearest 
  value to imputate missing values.
}
\examples{
 data(SurvLongData)
 # A truncated dataset to keep example run time brief
 exp <- nearValue(X = X[1:100,], Z = Z, tau = 1.0)

}
\references{
Cao H., Churpek M. M., Zeng D., Fine J. P. 
 (2015).
 Analysis of the proportional hazards model with sparse longitudinal covariates.
 Journal of the American Statistical Association, 110, 1187-1196.
}
\seealso{
\code{\link{fullKernel}}, \code{\link{halfKernel}}, \code{\link{lastValue}}
}

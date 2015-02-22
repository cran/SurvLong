\name{lastValue}
\alias{lastValue}
\title{
Last Value Carried Forward Method
}
\description{
A simple approach to evaluate the effects of longitudinal covariates 
on the occurrence of events when the time-dependent covariates are 
measured intermittently. 
Regression parameter are estimated using last value carried forward
imputation of missing values.
}
\usage{
lastValue(X, Z, tau, tol = 0.001, maxiter = 100 )
}
\arguments{
  \item{X    }{an object of class data.frame. 
               The structure of the data.frame must be 
               \{patient ID, event time, event indicator\}. 
               Patient IDs must be of class integer or be able to be coerced to 
               class integer without loss of information.
               Missing values must be indicated as NA.
                The event indicator is 1 if the event occurred; 0 if censored.}
  \item{Z    }{an object of class data.frame. 
               The structure of the data.frame must be 
               \{patient ID, time of measurement, measurement(s)\}. 
               Patient IDs must be of class integer or be able to be coerced to 
               class integer without loss of information.
               Missing values must be indicated as NA.}
  \item{tau  }{an object of class numeric. The desired time point.}
  \item{tol  }{An object of class numeric. The minimum change in the
               regression parameters deemed to indicate convergence of the
               Newton-Raphson method.}
  \item{maxiter }{ An object of class numeric. The maximum number of iterations
                  used to estimate regression parameters.}
}
\value{
A list is returned. 

  \item{betaHat }{The estimated model coefficients.}
  \item{stdErr  }{The standard error for each coefficient.}
  \item{zValue  }{The estimated z-value for each coefficient.}
  \item{pValue  }{The p-value for each coefficient.}

}
\references{
Cao, H., Churpek, M. M., Zeng, D., and Fine, J. P. (2015).
Journal of the American Statistical Association, in press.
}
\author{
Hongyuan Cao, Mathew M. Churpek, Donglin Zeng, Jason P. Fine, and Shannon T. Holloway
}
\seealso{
\code{\link{fullKernel}}, \code{\link{halfKernel}}, \code{\link{nearValue}}
}
\examples{
  data(SurvLongData)
  # A truncated dataset to keep example run time brief
  exp <- lastValue(X = X[1:200,], Z = Z, tau = 1.0)

}


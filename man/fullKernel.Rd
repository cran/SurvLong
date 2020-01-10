\name{fullKernel}
\alias{fullKernel}
\title{
Full Kernel Estimation with Forward and Backward Lagged Covariates
}
\description{
A kernel weighting scheme to evaluate the effects of longitudinal covariates 
on the occurrence of events when the time-dependent covariates are 
measured intermittently.
Regression parameter estimation uses full kernel imputation of missing 
values with both forward and backward lagged covariates.
}
\usage{
fullKernel(X, Z, tau, kType = "epan", bw = NULL, tol = 0.001, 
           maxiter = 100, verbose = TRUE)
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
  \item{kType}{An object of class character indicating the type of 
               smoothing kernel to use in the estimating equation. 
               Must be one of \{"epan", "uniform", "gauss"\}, where
               "epan" is the Epanechnikov kernel and "gauss" is the
               Gaussian kernel.}
  \item{bw   }{If provided, bw is an object of class numeric or a 
               numeric vector containing the bandwidths for which parameter
               estimates are to be obtained. 
               If NULL, an optimal bandwidth will be determined
               using an adaptive selection procedure. 
               The range of the bandwidth search space is taken
               to be 2*(Q3 - Q1)*n^{-0.7} to 2*(Q3 - Q1)*n^{-0.3},
               where Q3 is the 0.75 quantile and Q1 is the 0.25 quantile
               of the measurement times for the covariate
               and n is the effective number of patients, taken
               as the total number of patients that experienced an event.  }
  \item{tol  }{An object of class numeric. The minimum change in the
               regression parameters deemed to indicate convergence of the
               Newton-Raphson method.}
  \item{maxiter }{ An object of class numeric. The maximum number of iterations
                  used to estimate regression parameters.}
  \item{verbose }{ An object of class logical. TRUE results in screen prints.}
}
\value{
A list is returned. If bandwidths are provided, each element 
of the list is a matrix,
where the ith row corresponds to the ith bandwidth of input argument ``bw," 
and the columns correspond to the model parameters. If the bandwidth is
determined automatically, each element of the list is a named vector calculated
at the optimal bandwidth.

  \item{betaHat }{The estimated model coefficients.}
  \item{stdErr  }{The standard error for each coefficient.}
  \item{zValue  }{The estimated z-value for each coefficient.}
  \item{pValue  }{The p-value for each coefficient.}

If the bandwidth is determined automatically, three additional list
elements are returned:

  \item{optBW   }{The estimated optimal bandwidth.}
  \item{minMSE  }{The mean squared error at the optimal bandwidth.}
  \item{MSE     }{The vector of MSE for each bandwidth.}
}
\references{
Cao H., Churpek M. M., Zeng D., Fine J. P. 
(2015).
Analysis of the proportional hazards model with sparse longitudinal covariates.
Journal of the American Statistical Association, 110, 1187-1196.
}
\author{
Hongyuan Cao, Mathew M. Churpek, Donglin Zeng, Jason P. Fine, and Shannon T. Holloway
}
\seealso{
\code{\link{halfKernel}}, \code{\link{lastValue}}, \code{\link{nearValue}}
}
\examples{

  data(SurvLongData)

  exp <- fullKernel(X = X, Z = Z, tau = 1.0, bw = 0.015)

}


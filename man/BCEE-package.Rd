\name{BCEE-package}
\alias{BCEE-package}
\alias{BCEE}
\docType{package}
\title{
The Bayesian Causal Effect Estimation (BCEE) Algorithm}
\description{
Contains two functions that correspond to the approximate BCEE (A-BCEE) and the naive BCEE (N-BCEE) implementations of the BCEE algorithm described in Talbot et al. (2015).
}
\details{
\tabular{ll}{
Package: \tab BCEE\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-06-24\cr
License: \tab GPL (>=2) \cr
}
ABCEE(X, Y, U, omega), \cr
NBCEE(X, Y, U, omega).
}
\author{
Denis Talbot, Genevieve Lefebvre, Juli Atherton.

Maintainer: Denis Talbot \email{denis.talbot@fmed.ulaval.ca}
}
\references{
Talbot, D., Lefebvre, G., Atherton, J. (2015) \emph{The Bayesian causal effect estimation algorithm}, Journal of Causal Inference (in press).
}
\keyword{causal}
\keyword{confounding}
\keyword{model average}
\seealso{
 \code{\link{ABCEE}}, \code{\link{NBCEE}}.
}

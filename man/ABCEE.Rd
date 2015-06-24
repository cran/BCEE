\name{ABCEE}
\alias{ABCEE}
\title{
Approximate BCEE Implementation
}
\description{
A-BCEE implementation of the BCEE algorithm.
}
\usage{
ABCEE(X, Y, U, omega, forX = NA, niter = 5000, nburn = 500, nthin = 10,
 maxmodelY = NA, OR = 20)
}
\arguments{
  \item{X}{A vector of observed values for the continuous exposure.
}
  \item{Y}{A vector of observed values for the continuous outcome.
}
  \item{U}{A matrix of observed values for the \code{M} potential confounding covariates, where each column contains observed values for a potential confounding factor. A recommended implementation is to only consider pre-exposure covariates.
}
  \item{omega}{The value of the hyperparameter omega in the BCEE's outcome model prior distribution. A recommended implementation is to take \code{omega}\code{ = }\code{sqrt(n)*c}, where \code{n} is the sample size and \code{c} is a user-supplied constant value. Simulation studies suggest that values of \code{c} between 100 and 1000 yield good results. 
}
  \item{forX}{
A Boolean vector of size \code{M}, where the \code{m}th element indicates whether or not the \code{m}th potential confounding covariate should be considered in the exposure modeling step of the BCEE algorithm. The default for \code{forX} is \code{NA}, which indicates that all potential confounding covariates should be considered in the exposure modeling step. 
}
  \item{niter}{The number of post burn-in iterations in the Markov chain Monte Carlo model composition (MC^3) algorithm (Madigan et al. 1995), prior to applying thinning. The default is 5000.
}
  \item{nburn}{
The number of burn-in iterations (prior to applying thinning). The default is 500.
}
  \item{nthin}{
The thinning of the chain. The default is 10. 
}
  \item{maxmodelY}{
The maximum number of distinct outcome models that the algorithm can explore. Choosing a smaller value can shorten computing time. However, choosing a value that is too small will cause the algorithm to crash. The default is \code{NA}; the maximum number of outcome models that can be explored is then set to the minimum of \code{niter + nburn} and \code{2^M}.  
}
  \item{OR}{
A number specifying the maximum ratio for excluding models in Occam's window for the exposure modeling step (see the \code{bic.glm} help file, and Madigan & Raftery, 1994).  The default is 20.  
}
}

\details{
The \code{ABCEE} function first computes the exposure model's posterior distribution using the \code{bic.glm} function. The outcome model's posterior distribution is then computed using MC^3 (Madigan et al., 1995) as described in Talbot et al. (2015).

\code{ABCEE} assumes there are no missing values in the objects \code{X}, \code{Y} and \code{U}. The \code{na.omit} function which removes cases with missing data or an imputation package might be helpful.
}

\value{
\item{betas}{A vector containing the sampled values for the exposure effect.} 
\item{models.X}{A matrix giving the posterior distribution of the exposure model. Each row corresponds to an exposure model. Within each row, the first \code{M} elements are Booleans indicating the inclusion (1) or the exclusion (0) of each potential confounding factor. The last element gives the posterior probability of the exposure model.}
\item{models.Y}{A Boolean matrix identifying the sampled outcome models. Each row corresponds to a sampled outcome model. Within each row, the \code{m}th element equals 1 if and only if the \code{m}th potential confounding covariate is included in the sampled outcome model (and 0 otherwise).}
}

\references{
Madigan, D., York, J., Allard, D. (1995) \emph{Bayesian graphical models for discrete data}, International Statistical Review, 63, 215-232.

Madigan, D., Raftery, A. E. (1994) \emph{Model selection and accounting for model uncertainty in graphical models using Occam's window}, Journal of the American Statistical Association, 89 (428), 1535-1546.

Talbot, D., Lefebvre, G., Atherton, J. (2015) \emph{The Bayesian causal effect estimation algorithm}, Journal of Causal Inference (in press).
}
\author{
Denis Talbot, Genevieve Lefebvre, Juli Atherton.
}

\seealso{
\code{\link[BMA]{bic.glm}}, \code{\link{na.omit}}, \code{\link{NBCEE}}.
}
\examples{
# In this example, U1 and U2 are potential confounding covariates
# generated as independent N(0,1).
# X is generated as a function of both U1 and U2 with a N(0,1) error.
# Y is generated as a function of X and U1 with a N(0,1) error.
# Variable U1 is the only confounder.
# The causal effect of X on Y equals 1. 
# The exposure effect estimator (beta hat) in the outcome model  
# including U1 and U2 or including U1 only is unbiased.
# The sample size is n = 200.

# Generating the data
set.seed(418949); 
U1 = rnorm(200); 
U2 = rnorm(200);
X = 0.5*U1 + 1*U2 + rnorm(200);
Y = 1*X + 0.5*U1 + rnorm(200);

# Using ABCEE to estimate the causal exposure effect
n = 200;
omega.c = 500;
results = ABCEE(X,Y,cbind(U1,U2), omega = omega.c*sqrt(n),
 niter = 1000, nthin = 5, nburn = 20);

# The posterior mean of the exposure effect:
mean(results$betas);
# The posterior standard deviation of the exposure effect:
sd(results$betas);
# The posterior distribution of the exposure model:
results$models.X;
# The posterior probability of inclusion of each covariate:
colMeans(results$models.Y);
# The posterior distribution of the outcome model:
table(apply(results$models.Y, 1, paste0, collapse = ""));
}
\keyword{causal}
\keyword{confounding}
\keyword{model average}

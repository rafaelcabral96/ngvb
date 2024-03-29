% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ngvb.R
\name{ngvb}
\alias{ngvb}
\title{Approximate inference of latent non-Gaussian models (LnGMs)}
\usage{
ngvb(
  fit = NULL,
  manual.configs = NULL,
  selection = NULL,
  alpha.eta = rep(1, length(selection)),
  n.sampling = 1000,
  d.sampling = FALSE,
  method = "SCVI",
  fast = FALSE,
  verbose = TRUE,
  history = FALSE,
  V.init = NULL,
  eta.init = NULL,
  iter = 10,
  stop.rel.change = NULL
)
}
\arguments{
\item{fit}{The inla object that fits the LGM.}

\item{manual.configs}{Not necessary if \code{fit} is provided. List containing:
\itemize{
   \item \code{inla.fit.V} Function that receives as input the mixing variables \eqn{\mathbf{V}} and
   the output is an inla object that fits the LGM \eqn{(\mathbf{x},\boldsymbol{\theta} | \mathbf{V}, \mathbf{y})}.
   \item \code{Dfunc} Function that receives as input the hyperparameter vector \eqn{\boldsymbol{\theta}}
    in internal scale and the output is \eqn{\mathbf{D}(\boldsymbol{\theta})}, where \eqn{\mathbf{D}} is the dependency matrix that
    specifies the non-Gaussian latent field. If there is more than one latent component
    to be extended to non-Gaussianity, this should be a list of functions \eqn{\mathbf{D}_i(\boldsymbol{\theta})},
    where \eqn{\mathbf{D}_i(\boldsymbol{\theta})} is the dependency matrix that specifies component i.
    \item \code{h} Predefined constant \eqn{\mathbf{h}} vector that contains the distance between locations
    or area of the basis functions. For models defined in discrete space this should be a
    vector of ones. If there is more than one latent component to be extended to non-Gaussianity,
    this should be a list of vectors \eqn{\mathbf{h_i}}
    where \eqn{\mathbf{h_i}} is the predefined constant vector of component i.
 }}

\item{selection}{List which specifies which model components of the LGM are to be extended
to non-Gaussianity. Same syntax as the argument \code{selection} of
the function \code{inla.posterior.sample}.}

\item{alpha.eta}{Numeric. Rate parameter of the exponential prior of the non-Gaussianity parameter \eqn{\alpha}.
Should be a vector with the same dimension as the number of model components to extend.}

\item{n.sampling}{Numeric. Number of samples uses in several sampling task throughout the SVI or SCVI algorithm.}

\item{d.sampling}{Logical. If \code{TRUE} the expectations \eqn{d_i = E([\mathbf{D}\mathbf{x}]_i^2)} are computed
by sampling (slower). If \code{FALSE} (default) uses the mixture Gaussian approximation
of \eqn{\mathbf{x}} obtained in INLA.}

\item{method}{Character. Should be \code{"SCVI"}, \code{"SVI"}, or \code{"Gibbs"} if the algorithm
should be the structural and collapsed variational inference algorithm (faster) or the structural
variational inference algorithm (slower), or a Gibbs sampler (even slower), respectively.}

\item{fast}{Logical. If \code{TRUE} then INLA will be run with the empirical Bayes mode and
several control variables will be turned off to improve the speed of the algorithm.}

\item{verbose}{Logical. If \code{TRUE} print the posterior mean of \eqn{\eta} and plot the posterior
mean of \eqn{\mathbf{V}} at each iteration.}

\item{history}{Logical. If \code{TRUE} save LGM inla object of each iteration. If \code{FALSE} only
save the LGM inla object of the last iteration.}

\item{V.init}{List of numeric vectors. Initial values for \eqn{\mathbf{V}}. By default \eqn{\mathbf{V}=\mathbf{h}}.}

\item{eta.init}{Numeric vector. Initial values of \eqn{\eta}. By default eta.init = 0.5.}

\item{iter}{Integer. Maximum number of iterations.}

\item{stop.rel.change}{Numeric. Stop the algorithm when the relative change in the posterior mean
of \eqn{\eta} is smaller than \code{stop.rel.change}.}
}
\value{
An S4 object containing the outputs. Slots are accessed with @. These are:
\itemize{
   \item \code{history} List containing summaries and inla objects (if \code{history = TRUE}) of each iteration (except the last one).
   \item \code{LGM} inla object containing summaries and marginals of the latent field \eqn{\mathbf{x}}
   and hyperparameters \eqn{\boldsymbol{\theta}} of the last iteration.
   \item \code{summary.mixing} Data frame containing summaries of the mixing variables \eqn{\mathbf{V}}.
   \item \code{summary.ng} Data frame containing summaries of the non-Gaussianiy parameters \eqn{\eta}.
   \item \code{configs} List containing the LnGM model specifications
 }
}
\description{
Finds the posterior distribution of the hyperparameters \eqn{\boldsymbol{\theta}}, latent field \eqn{\mathbf{x}},
mixing variables \eqn{\mathbf{V}}, and non-Gaussianity parameters \eqn{\eta} of latent non-Gaussian models (LnGMs)
using INLA and variational Bayes approximations. These LnGMs are flexible extensions
of latent Gaussian models (LGMs). The LnGMs are specified either through
an INLA object (\code{fit}) that fits an LGM or through a list of configurations (\code{manual.configs}).
Run `devtools::build_vignettes("ngvb")` and `vignette("ngvb")` to see several use cases.
}
\examples{
\dontrun{
 #Here we fit an RW1 latent process to the jumpts time series
 plot(jumpts)

 #Fit LGM with INLA
 LGM     <- inla(y ~ -1 + f(x,  model = "rw1"),
                 data = jumpts)

 #Check adequacy of latent Gaussianity assumption
 check.list <- ng.check(fit = LGM)

 #Fit LnGM with ngvb
 LnGM <- ngvb(fit = LGM)

 #Available methods
 summary(LnGM)
 print(LnGM)
 plot(LnGM)
 fitted(LnGM)
 samples <- simulate(LnGM)
 }
}

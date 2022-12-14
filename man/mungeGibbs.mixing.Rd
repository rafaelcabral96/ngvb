% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ngvb.list.R
\name{mungeGibbs.mixing,ngvb.list-method}
\alias{mungeGibbs.mixing,ngvb.list-method}
\title{Process \code{ngvb.list} when \code{method = "Gibbs"}. Produces a matrix with the
samples of the mixing variables \eqn{\mathbf{V}}.}
\usage{
\S4method{mungeGibbs.mixing}{ngvb.list}(object)
}
\arguments{
\item{object}{An ngvb.list object (output of \code{ngvb} function)}
}
\description{
Process \code{ngvb.list} when \code{method = "Gibbs"}. Produces a matrix with the
samples of the mixing variables \eqn{\mathbf{V}}.
}
\examples{
\dontrun{
 #Here we fit an RW1 latent process to the jumpts time series
 plot(jumpts)

 #Fit LGM with INLA
 LGM     <- inla(y ~ -1 + f(x,  model = "rw1"),
                 data = jumpts)

 #Run 10 iterations from the Gibbs sampler
 LnGM.Gibbs <- ngvb(fit = LGM, selection = list(x=1:100),
                    method = "Gibbs", iter=10)

 Gibbs.V    <- mungeGibbs.mixing(LnGM.Gibbs)
 Gibbs.eta  <- mungeGibbs.ng(LnGM.Gibbs)
 }
}

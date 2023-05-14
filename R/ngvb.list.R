#' This class stores the output of the ngvb function.
#'
#' @slot history List containing summaries of the model at each iteration.
#' @slot LGM  inla object containing summaries and marginals of the latent field \eqn{\mathbf{x}}
#'    and hyperparameters \eqn{\boldsymbol{\theta}} of the last iteration.
#' @slot summary.mixing  Data frame containing summaries of the mixing variables \eqn{\mathbf{V}} (last iteration).
#' @slot summary.ng Data frame containing summaries of the non-Gaussianiy parameters \eqn{\eta} (last iteration).
#' @slot configs List containing the LnGM model specifications.
ngvb.list <- setClass("ngvb.list",
                      slots = c(history = "list",
                                LGM     = "inla",
                                summary.mixing = "data.frame",
                                summary.ng = "data.frame",
                                configs = "list"))


setGeneric("addfit",
           function(fits, fit, history) StandardGeneric("addfit"))

setMethod("addfit", signature("ngvb.list","inla","logical"),
          function(fits, fit, history){

            n <- length(fits@history)

            #if not the first iteration then add old elements to history
            if(!is.null(fit$ngvb$d)){
              fits@history[[n+1]] <- list(summary.mixing = fits@summary.mixing, summary.ng = fits@summary.ng)
              if(history) fits@history[[n+1]] <- append(fits@history[[n+1]], LGM = fits@LGM)
            }

            fits@LGM            <- fit
            fits@summary.mixing <- ngvb.summary.mixing(fit, fits@configs)
            fits@summary.ng     <- ngvb.summary.ng(fit, fits@configs)

            return(fits)
          }
)


setGeneric("addconfigs",
           function(fits, list) StandardGeneric("addconfigs"))

setMethod("addconfigs", signature("ngvb.list","list"),
          function(fits, list){

            fits@configs    <- list

            return(fits)
          }
)

setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

#' Plots several summaries of an \code{ngvb.list} object.
#'
#' @param x An ngvb.list object (output of \code{ngvb} function)
#' @param y Not used.
#' @param ... Extra arguments to be used in \code{plot(x@LGM, ...)} where \code{x@LGM} is an inla object.
#' @examples
#' \dontrun{
#'  #Here we fit an RW1 latent process to the jumpts time series
#'  plot(jumpts)
#'
#'  #Fit LGM with INLA
#'  LGM     <- inla(y ~ -1 + f(x,  model = "rw1"),
#'                  data = jumpts)
#'
#'  #Fit LnGM with ngvb
#'  LnGM <- ngvb(fit = LGM, selection = list(x=1:100))
#'
#'  #Available methods
#'  summary(LnGM)
#'  print(LnGM)
#'  plot(LnGM)
#'  fitted(LnGM)
#'  samples <- simulate(LnGM)
#'  }
#' @export
#' @rdname plot
setMethod("plot",
          c(x="ngvb.list", y="missing"),
          function(x, y, ...){

            #Inla plots
            #plot(x@LGM, ...)


            #Mixing variables
            groups <- unique(x@summary.mixing$component)
            par(mfrow=c(1,1))
            for(group in groups){
              sel <- x@summary.mixing$component == group
              plot(x@summary.mixing$mean[sel], xlab = "index", ylab = "mean", main = paste0("Mixing variables for ", group))
            }

            #Evolution of eta
            prog <- append(lapply(x@history, function(y) y$summary.ng$mean),list(x@summary.ng$mean))
            prog <- as.data.frame(do.call(rbind, prog))
            colnames(prog) <- paste0("mean for ", groups)

            plot.ts(prog, main = "Evolution of the non-Gaussianity parameters", xlab = "iteration")

            #Marginal of eta

          })

#' Produces summaries of the \code{ngvb.list} object.
#'
#' @param object An ngvb.list object (output of \code{ngvb} function)
#' @examples
#' \dontrun{
#'  #Here we fit an RW1 latent process to the jumpts time series
#'  plot(jumpts)
#'
#'  #Fit LGM with INLA
#'  LGM     <- inla(y ~ -1 + f(x,  model = "rw1"),
#'                  data = jumpts)
#'
#'  #Fit LnGM with ngvb
#'  LnGM <- ngvb(fit = LGM, selection = list(x=1:100))
#'
#'  #Available methods
#'  summary(LnGM)
#'  print(LnGM)
#'  plot(LnGM)
#'  fitted(LnGM)
#'  samples <- simulate(LnGM)
#'  }
#' @export
#' @rdname summary
setMethod("summary", "ngvb.list",
          function(object){

            sum.LGM <- summary(object@LGM)
            sum.LGM$cpu.used       <- NA
            sum.LGM$total.time     <- object@LGM$ngvb$time
            sum.LGM$summary.ng     <- object@summary.ng
            sum.LGM$summary.mixing <- object@summary.mixing

            print(object@summary.ng)

            return(sum.LGM)
          })

#' Prints outputs of the \code{ngvb.list} object.
#'
#' @param x An ngvb.list object (output of \code{ngvb} function)
#' @examples
#' \dontrun{
#'  #Here we fit an RW1 latent process to the jumpts time series
#'  plot(jumpts)
#'
#'  #Fit LGM with INLA
#'  LGM     <- inla(y ~ -1 + f(x,  model = "rw1"),
#'                  data = jumpts)
#'
#'  #Fit LnGM with ngvb
#'  LnGM <- ngvb(fit = LGM, selection = list(x=1:100))
#'
#'  #Available methods
#'  summary(LnGM)
#'  print(LnGM)
#'  plot(LnGM)
#'  fitted(LnGM)
#'  samples <- simulate(LnGM)
#'  }
#' @export
#' @rdname print
setMethod(f = "print", "ngvb.list",
          function(x){
            print(x@LGM)
            cat("\n \nngvb Configurations: \n")
            print(x@configs)

          })

#' Summary of the fitted values. Also returns the marginals of the fitted values
#' if \code{compute=TRUE} in \code{control.predictor} of the inla object \code{fit}.
#'
#' @param object An ngvb.list object (output of \code{ngvb} function)
#' @examples
#' \dontrun{
#'  #Here we fit an RW1 latent process to the jumpts time series
#'  plot(jumpts)
#'
#'  #Fit LGM with INLA
#'  LGM     <- inla(y ~ -1 + f(x,  model = "rw1"),
#'                  data = jumpts)
#'
#'  #Fit LnGM with ngvb
#'  LnGM <- ngvb(fit = LGM, selection = list(x=1:100))
#'
#'  #Available methods
#'  summary(LnGM)
#'  print(LnGM)
#'  plot(LnGM)
#'  fitted(LnGM)
#'  samples <- simulate(LnGM)
#'  }
#' @export
#' @rdname fitted
setMethod(f = "fitted", "ngvb.list",
          function(object){
            fitted <- list(summary  = object@LGM$summary.fitted.values,
                           marginal = object@LGM$marginals.fitted.values)
            return(fitted)
          })



#' Obtain samples from the fitted LnGM model parameters.
#' @param object An ngvb.list object (output of \code{ngvb} function).
#' @param n.sampling Integer. Number of samples.
#' @param components Vector containing \code{c("LGM","V","ng","hyperpar")}.
#' If \code{components} contains \code{"LGM"} then it generates samples of \eqn{(\mathbf{x},\boldsymbol{\theta})} using
#' \code{inla.posterior.sample}. If it contains \code{"hyperpar"} then generate samples of
#' \eqn{\boldsymbol{\theta}} using \code{inla.hypearpar.sample}. It it contains \code{"V"} or \code{"ng"} then
#' generate samples from the mixing variables \eqn{\mathbf{V}} and \eqn{\eta} for each model component,
#' respectively.
#' @param improved.tail Logical. If \code{TRUE} generate leptokurtic samples of the the latent field
#' \eqn{\mathbf{x}}. It first generates samples of \eqn{(\mathbf{V},\boldsymbol{\theta})} and then
#' it generates samples of  \eqn{\mathbf{x} | \mathbf{V} , \boldsymbol{\theta}, \mathbf{y}}
#' (which is an LGM) by fitting an INLA model each time and generating n = \code{augmentation} samples.
#'  Slow. Reduce \code{n.sampling} for speed and increase \code{augmentation} to obtain more samples.
#' @param augmentation Integer. If \code{improved.tail = TRUE}, then for each sample of \eqn{(\mathbf{V},\boldsymbol{\theta})} generate
#' n = \code{augmentation} samples of \eqn{\mathbf{x} | \mathbf{V} , \boldsymbol{\theta} , \mathbf{y}}.
#' @examples
#' \dontrun{
#'  #Here we fit an RW1 latent process to the jumpts time series
#'  plot(jumpts)
#'
#'  #Fit LGM with INLA
#'  LGM     <- inla(y ~ -1 + f(x,  model = "rw1"),
#'                  data = jumpts)
#'
#'  #Fit LnGM with ngvb
#'  LnGM <- ngvb(fit = LGM, selection = list(x=1:100))
#'
#'  #Available methods
#'  summary(LnGM)
#'  print(LnGM)
#'  plot(LnGM)
#'  fitted(LnGM)
#'  samples <- simulate(LnGM)
#'  }
#' @export
#' @rdname simulate
setMethod(f = "simulate", "ngvb.list",
          function(object, n.sampling = object@configs$n.sampling,
                   components = c("LGM","V","ng"), improved.tail = FALSE,
                   augmentation = 10){

            ncomp       <- object@configs$ncomp
            method      <- object@configs$method
            N           <- object@configs$N

            hyper.sample <- list()
            LGM.sample   <- list()
            eta.sample   <- list()
            V.sample     <- list()

            if(improved.tail) components <- c(components, "V", "hyperpar")

            if("hyperpar" %in% components){
              hyper.sample <- inla.hyperpar.sample(n.sampling, object@LGM, improve.marginals = TRUE, inter=TRUE)

            }

            if("V" %in% components){
              if(method == "SCVI"){
                if(n.sampling > nrow(object@LGM$ngvb$SCVI.samples[[1]]) ){

                  d    <- object@LGM$ngvb$d
                  h    <- object@configs$h
                  theta_eta <- object@configs$theta_eta
                  N <- object@configs$N

                  for(k in 1:ncomp){
                    dk <- d[[k]]
                    hk <- h[[k]]
                    theta_etak <- theta_eta[k]
                    Nk    <- N[k]

                    eta.prior <- eta.prior.f(dk,hk,theta_etak,Nk)
                    eta.s <- sampler.inverseCDF(eta.prior, supp.min = 0, supp.max = 1000, supp.points = 10000, n.samples = n.sampling)

                    ntotal <- n.sampling*length(hk)
                    p      <- rep(-1,ntotal)
                    a      <- rep(1/eta.s, each=length(hk))
                    b      <- c(t(outer(eta.s, hk, FUN = function(x,y) y^2/x))) + rep(dk,length(eta.s))
                    Vm     <- matrix(1/ngme::rGIG(p,a,b), nrow = n.sampling, ncol = length(hk), byrow = TRUE)

                    V.sample[[k]]   <- 1/Vm
                    eta.sample[[k]] <- eta.s

                  }

                }
                else{
                  V.sample <- object@LGM$ngvb$SCVI.samples
                }
              }
              else if(method == "SVI"){

                h <- object@configs$h
                N <- object@configs$N

                comps <- names(object@configs$selection)
                for(k in 1:ncomp){
                  hk    <- h[[k]]
                  Nk    <- N[k]
                  summary.mixing.k = object@summary.mixing[object@summary.mixing$component == comps[k],]
                  p      <- rep(object@summary.mixing$GIG.p, each = n.sampling)
                  a      <- rep(object@summary.mixing$GIG.a, each = n.sampling)
                  b      <- rep(object@summary.mixing$GIG.b, each = n.sampling)
                  Vm     <- matrix(ngme::rGIG(p,a,b), nrow = n.sampling, ncol = length(hk), byrow = FALSE)
                  V.sample[[k]]   <- Vm
                }
              }

              names(V.sample) <- names(object@configs$selection)
            }

            if("ng" %in% components){
              if(method == "SCVI"){
                if(n.sampling >length(object@LGM$ngvb$eta.samples[[1]])){
                  if("V" %in% components){
                    #do nothing
                  }
                  else{
                    d    <- object@LGM$ngvb$d
                    h    <- object@configs$h
                    theta_eta <- object@configs$theta_eta
                    N <- object@configs$N

                    for(k in 1:ncomp){
                      dk <- d[[k]]
                      hk <- h[[k]]
                      theta_etak <- theta_eta[k]
                      Nk    <- N[k]

                      eta.prior <- eta.prior.f(dk,hk,theta_etak,Nk)
                      eta.sample[[k]] <- sampler.inverseCDF(eta.prior, supp.min = 0, supp.max = 1000, supp.points = 10000, n.samples = n.sampling)
                    }
                  }
                }else{
                  eta.sample <- object@LGM$ngvb$eta.samples
                }
              }

              else if(method == "SVI"){
                for(k in 1:ncomp){
                  summary.mixing.k = object@summary.ng[k,]
                  p      <- rep(summary.mixing.k$GIG.p, each = n.sampling)
                  a      <- rep(summary.mixing.k$GIG.a, each = n.sampling)
                  b      <- rep(summary.mixing.k$GIG.b, each = n.sampling)
                  eta.sample[[k]]     <- ngme::rGIG(p,a,b)
                }
              }

              names(eta.sample) <- names(object@configs$selection)

            }

            if("LGM" %in% components){
              selection   <- object@configs$selection
              fit        <- object@LGM

              if(improved.tail){

                inla.fit.V <- object@configs$inla.fit.V

                opts <- list()
                pb <- txtProgressBar(min=1, max=n.sampling, style=3)
                progress <- function(n) setTxtProgressBar(pb, n)
                opts$progress <- progress


                #For a given value of theta and V, fit INLA object and find post distribution of x|y,theta,V
                samples <- list()
                cat(paste("0Oo----------- Generating samples from latent field x  -----------oO0\n", sep = ""))
                LGM.sample <- foreach(i = 1:n.sampling, .combine = 'c',
                                      .packages = c("INLA"),
                                      .export = c("inla.fit.V","fit","hyper.sample","V.sample","selection"),
                                      .options.snow=opts,
                                      .errorhandling = 'pass') %do% {

                                        theta    <- hyper.sample[i,]
                                        V        <- lapply(V.sample,'[',i,)
                                        inla.fit <- inla.fit.V(V,theta)
                                        inla.posterior.sample(n = augmentation, inla.fit, selection = selection,
                                                              intern = FALSE, use.improved.mean = TRUE, skew.corr = TRUE,
                                                              add.names = TRUE, seed = 0L, num.threads = 0,
                                                              verbose = FALSE)
                                      }


              }else{

                LGM.sample <- inla.posterior.sample(n = n.sampling, fit, selection = selection,
                                                    intern = TRUE, use.improved.mean = TRUE, skew.corr = TRUE,
                                                    add.names = TRUE, seed = 0L, num.threads = 0,
                                                    verbose = FALSE)
              }

            }


            return(list(hyperpar = hyper.sample, LGM = LGM.sample, V = V.sample, ng = eta.sample))


          })



setGeneric("mungeGibbs.mixing",
           function(object) StandardGeneric("mungeGibbs.mixing"))

#' Process \code{ngvb.list} when \code{method = "Gibbs"}. Produces a matrix with the
#' samples of the mixing variables \eqn{\mathbf{V}}.
#'
#' @param object An ngvb.list object (output of \code{ngvb} function)
#' @examples
#' \dontrun{
#'  #Here we fit an RW1 latent process to the jumpts time series
#'  plot(jumpts)
#'
#'  #Fit LGM with INLA
#'  LGM     <- inla(y ~ -1 + f(x,  model = "rw1"),
#'                  data = jumpts)
#'
#'  #Run 10 iterations from the Gibbs sampler
#'  LnGM.Gibbs <- ngvb(fit = LGM, selection = list(x=1:100),
#'                     method = "Gibbs", iter=10)
#'
#'  Gibbs.V    <- mungeGibbs.mixing(LnGM.Gibbs)
#'  Gibbs.eta  <- mungeGibbs.ng(LnGM.Gibbs)
#'  }
#' @export
#' @rdname mungeGibbs.mixing
setMethod("mungeGibbs.mixing", "ngvb.list",
          function(object){

            names <- names(object@configs$selection)
            n     <- length(object@history)
            ncomp <- object@configs$ncomp
            N     <- object@configs$N

            Vlist <- list()
            for(k in 1:ncomp){
              Vmatrix <- matrix(NA, nrow = n+1, ncol = N[k])
              for(i in 1:n){
                sum <- object@history[[i]]$summary.mixing
                V <- sum[sum$component == names[k],"V"]
                Vmatrix[i,] <-  V
              }
              sum <- object@summary.mixing
              V <- sum[sum$component == names[k], "V"]
              Vmatrix[i+1,] <-  V

              Vlist[[k]] <- Vmatrix
            }
            names(Vlist) <- names

            return(Vlist)
          })

setGeneric("mungeGibbs.ng",
           function(object) StandardGeneric("mungeGibbs.ng"))

#' Process \code{ngvb.list} when \code{method = "Gibbs"}. Produces a matrix with the
#' samples of the non-Gaussianity parameter \eqn{\eta}.
#'
#' @param object An ngvb.list object (output of \code{ngvb} function)
#' @examples
#' \dontrun{
#'  #Here we fit an RW1 latent process to the jumpts time series
#'  plot(jumpts)
#'
#'  #Fit LGM with INLA
#'  LGM     <- inla(y ~ -1 + f(x,  model = "rw1"),
#'                  data = jumpts)
#'
#'  #Run 10 iterations from the Gibbs sampler
#'  LnGM.Gibbs <- ngvb(fit = LGM, selection = list(x=1:100),
#'                     method = "Gibbs", iter=10)
#'
#'  Gibbs.V    <- mungeGibbs.mixing(LnGM.Gibbs)
#'  Gibbs.eta  <- mungeGibbs.ng(LnGM.Gibbs)
#'  }
#' @export
#' @rdname mungeGibbs.ng
setMethod("mungeGibbs.ng", "ngvb.list",
          function(object){

            names <- names(object@configs$selection)
            n     <- length(object@history)
            ncomp <- object@configs$ncomp
            N     <- object@configs$N

            nglist <- list()
            for(k in 1:ncomp){
              ng <- c()
              for(i in 1:n){
                eta <- object@history[[i]]$summary.ng[k,1]
                ng <- c(ng,eta)
              }
              ng <- c(ng, object@summary.ng[k,1])
              nglist[[k]] <- ng
            }
            names(nglist) <- names
            return(nglist)
          })

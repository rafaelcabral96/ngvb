#' @title Approximate inference of latent non-Gaussian models (LnGMs)
#'
#' @description Finds the posterior distribution of the hyperparameters \eqn{\boldsymbol{\theta}}, latent field \eqn{\mathbf{x}},
#' mixing variables \eqn{\mathbf{V}}, and non-Gaussianity parameters \eqn{\eta} of latent non-Gaussian models (LnGMs)
#' using INLA and variational Bayes approximations. These LnGMs are flexible extensions
#' of latent Gaussian models (LGMs). The LnGMs are specified either through
#' an INLA object (\code{fit}) that fits an LGM or through a list of configurations (\code{manual.configs}).
#' Run `devtools::build_vignettes("ngvb")` and `vignette("ngvb")` to see several use cases.
#'
#'
#' @param fit The inla object that fits the LGM.
#' @param manual.configs Not necessary if \code{fit} is provided. List containing:
#' \itemize{
#'    \item \code{inla.fit.V} Function that receives as input the mixing variables \eqn{\mathbf{V}} and
#'    the output is an inla object that fits the LGM \eqn{(\mathbf{x},\boldsymbol{\theta} | \mathbf{V}, \mathbf{y})}.
#'    \item \code{Dfunc} Function that receives as input the hyperparameter vector \eqn{\boldsymbol{\theta}}
#'     in internal scale and the output is \eqn{\mathbf{D}(\boldsymbol{\theta})}, where \eqn{\mathbf{D}} is the dependency matrix that
#'     specifies the non-Gaussian latent field. If there is more than one latent component
#'     to be extended to non-Gaussianity, this should be a list of functions \eqn{\mathbf{D}_i(\boldsymbol{\theta})},
#'     where \eqn{\mathbf{D}_i(\boldsymbol{\theta})} is the dependency matrix that specifies component i.
#'     \item \code{h} Predefined constant \eqn{\mathbf{h}} vector that contains the distance between locations
#'     or area of the basis functions. For models defined in discrete space this should be a
#'     vector of ones. If there is more than one latent component to be extended to non-Gaussianity,
#'     this should be a list of vectors \eqn{\mathbf{h_i}}
#'     where \eqn{\mathbf{h_i}} is the predefined constant vector of component i.
#'  }
#' @param selection List which specifies which model components of the LGM are to be extended
#' to non-Gaussianity. Same syntax as the argument \code{selection} of
#' the function \code{inla.posterior.sample}.
#' @param alpha.eta Numeric. Rate parameter of the exponential prior of the non-Gaussianity parameter \eqn{\alpha}.
#' Should be a vector with the same dimension as the number of model components to extend.
#' @param n.sampling Numeric. Number of samples uses in several sampling task throughout the SVI or SCVI algorithm.
#' @param d.sampling Logical. If \code{TRUE} the expectations \eqn{d_i = E([\mathbf{D}\mathbf{x}]_i^2)} are computed
#' by sampling (slower). If \code{FALSE} (default) uses the mixture Gaussian approximation
#' of \eqn{\mathbf{x}} obtained in INLA.
#' @param method Character. Should be \code{"SCVI"}, \code{"SVI"}, or \code{"Gibbs"} if the algorithm
#' should be the structural and collapsed variational inference algorithm (faster) or the structural
#' variational inference algorithm (slower), or a Gibbs sampler (even slower), respectively.
#' @param fast Logical. If \code{TRUE} then INLA will be run with the empirical Bayes mode and
#' several control variables will be turned off to improve the speed of the algorithm.
#' @param verbose Logical. If \code{TRUE} print the posterior mean of \eqn{\eta} and plot the posterior
#'  mean of \eqn{\mathbf{V}} at each iteration.
#' @param history Logical. If \code{TRUE} save LGM inla object of each iteration. If \code{FALSE} only
#' save the LGM inla object of the last iteration.
#' @param V.init List of numeric vectors. Initial values for \eqn{\mathbf{V}}. By default \eqn{\mathbf{V}=\mathbf{h}}.
#' @param eta.init Numeric vector. Initial values of \eqn{\eta}. By default eta.init = 0.5.
#' @param iter Integer. Maximum number of iterations.
#' @param stop.rel.change Numeric. Stop the algorithm when the relative change in the posterior mean
#' of \eqn{\eta} is smaller than \code{stop.rel.change}.
#' @return An S4 object containing the outputs. Slots are accessed with @@. These are:
#' \itemize{
#'    \item \code{history} List containing summaries and inla objects (if \code{history = TRUE}) of each iteration (except the last one).
#'    \item \code{LGM} inla object containing summaries and marginals of the latent field \eqn{\mathbf{x}}
#'    and hyperparameters \eqn{\boldsymbol{\theta}} of the last iteration.
#'    \item \code{summary.mixing} Data frame containing summaries of the mixing variables \eqn{\mathbf{V}}.
#'    \item \code{summary.ng} Data frame containing summaries of the non-Gaussianiy parameters \eqn{\eta}.
#'    \item \code{configs} List containing the LnGM model specifications
#'  }
#' @examples
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
#' @export
#'
ngvb <- function(fit = NULL, manual.configs = NULL,
                 selection, alpha.eta = rep(1,length(selection)),
                 n.sampling = 1000, d.sampling = FALSE,
                 method = "SCVI", fast = FALSE,
                 verbose = TRUE, history = FALSE,
                 V.init = NULL, eta.init = NULL,
                 iter = 10, stop.rel.change = NULL){


  start_time <- Sys.time()
  fits       <- ngvb.list()

  if(is.null(stop.rel.change)){
    if(method == "SCVI") stop.rel.change = 0.01
    if(method == "SVI") stop.rel.change = 0.001
  }

  ##-------------------------------------
  ##    HOW MANY COMPONENTS TO ROBUSTIFY?
  ##-------------------------------------
  ncomp = length(selection)

  ##-----------------------------------------------------
  ##    Getting Dfunc, h, and, inla.fit.V from fit object
  ##-----------------------------------------------------
  if(is.null(manual.configs)){
    fit                   <- change.args.controls(fit, fast)

    models                <- c("iid","rw1", "rw2", "ar1")

    ##-------------------------------------
    ##    CHECK SELECTION AND FIND MODEL TYPE
    ##-------------------------------------
    indexname      <- names(selection)
    #position of indexname in attr "term.label"
    pos            <- apply(as.array(indexname), MARGIN = 1, FUN = function(x) which(x == intersect(all.vars(fit$.args$formula),names(fit$.args$data))) - 1)
    term_string    <- attr(terms(fit$.args$formula), which = "term.label")[pos]
    #model          <- models[sapply(models, grepl, term_string)]
    model          <- apply(as.array(term_string), MARGIN =1, FUN = function(x) models[sapply(models, grepl, x)])
    #N = tail(selection[[indexname]], n=1) #dimension of latent field
    N              <- apply(as.array(indexname), MARGIN = 1, FUN = function(x) tail(selection[[x]], n=1) )
    if(length(model) < ncomp) stop("Model not implemented")

    ##-------------------------------------
    ##     Which hyperparameters are relevant?
    ##-------------------------------------
    hyperparsel       <- t(apply(as.array(indexname), MARGIN = 1,
                                 FUN = function(x) sapply(paste("\\b","for ", x, "\\b", sep=""), grepl, rownames(fit$summary.hyperpar))))

    ##-------------------------------------
    ##     Construct D.config - list with all configurations needed to compute discrepency measure d
    ##-------------------------------------
    D.config <- list()
    for(i in 1:ncomp){
      D.config[[i]] <-  find.D.config(fit, model[i], selection[i], indexname[i], hyperparsel[i,], n.sampling, N[i])
    }
    ##-------------------------------------
    ##    What is h?
    ##-------------------------------------
    h <- list()
    for(i in 1:ncomp){
      h[[i]]          <-  find.h(fit, D.config[[i]])
      D.config[[i]]$h <- h[[i]]}

    ##-------------------------------------
    ##    FIT INLA MODEL WITH FIXED V
    ##-------------------------------------
    inla.fit.V <- find.inla.fit.V(fit, N, model, term_string, ncomp, D.config)
  }

  else{
    ##---------------------------------------------------------
    ##    Getting Dfunc, h, and, inla.fit.V from manual.configs
    ##---------------------------------------------------------
    inla.fit.V   <- manual.configs$inla.fit.V
    indexname    <- names(selection)

    #manual.configs = list(inla.fit.V=3, Dfuncs = list(Dfunc=1,Dfunc=6), h =  list(h=3,h-4))
    Dfunc        <- manual.configs$Dfunc
    h            <- manual.configs$h

    N        <- list()
    D.config <- list()
    for(i in 1:ncomp){
      N[[i]]           <- length(selection[[i]])
      D.config[[i]]    <- list(indexname = names(selection)[i], Dfunc = Dfunc[[i]], h = h[[i]], n.sampling = n.sampling, selection = selection)}
  }

  ##-------------------------------------
  ##    First iteration already done
  ##-------------------------------------

  if(is.null(V.init)){
    V = h
  } else{V = V.init}

  fit             <- inla.fit.V(V)

  if(is.null(eta.init)){
    eta       <- rep(0.5, ncomp)
  }else{ eta  <- eta.init}



  fit$ngvb$V      <- V
  fit$ngvb$eta    <- eta
  fits <- addconfigs(fits, list(h = h, alpha.eta = alpha.eta, method = method,
                                N = N, ncomp = ncomp, selection = selection,
                                n.sampling = n.sampling,
                                inla.fit.V = inla.fit.V))

  fits <- addfit(fits, fit, history)

  cat(paste("\n0Oo----------- Iteration ", 1 ," -----------oO0\n", sep = ""))
  if(verbose){
    x11(width = ncomp*6, height = 6)
    par(mfrow=c(1,ncomp))
    for(k in 1:ncomp){
      plot(V[[k]]/h[[k]], main = paste("weights for", indexname[k]), ylab = "weights", xlab = "index")
    }

    cat("Initial value of eta: ")
    cat(round(eta,3))
    cat("\n")
  }


  if(iter>1) for(i in 2:iter){

    ##-------------------------------
    ##     Compute relevant statistic
    ##-------------------------------
    if(method == "SVI" || method == "SCVI"){
      if(d.sampling){d <- suppressWarnings(compute.d.sampling(fit, D.config)) }
      else{          d <- compute.d.configs(fit, D.config) }
    }


    if (method == "Gibbs"){

      samples <- inla.posterior.sample(n = 1, fit, selection = selection,
                                       intern = TRUE, use.improved.mean = TRUE, skew.corr = TRUE,
                                       add.names = TRUE, seed = 0L, num.threads = 0,
                                       verbose = FALSE)

      lenghts <- unlist(lapply(selection,length))
      start = 1
      end = lenghts[1]
      d <- list()
      for(k in 1:ncomp){

        alpha.etak <- alpha.eta[k]
        hk     <- h[[k]]
        etak   <- eta[k]
        Nk     <- N[[k]]

        Vk <- numeric(lenghts[k])

        x        = c(samples[[1]]$latent)[start:end]
        hyperpar = samples[[1]]$hyperpar
        Dx       = (D.config[[k]]$Dfunc(hyperpar)%*%x)[,1]

        d[[k]] <- Dx

        for(j in 1:length(hk)){
          #Vk[j] <- ngme::rGIG(p = -1, a = 1/etak, b = hk[j]^2/etak + Dx[j]^2 )
          Vk[j] <- GIGrvg::rgig(n = 1,
                                lambda = -1,
                                psi = 1/eta[k],
                                chi = hk[j]^2/eta[k] + Dx[j]^2)
        }

        V[[k]]       <- Vk
        eta[k]       <- GIGrvg::rgig(n = 1,
                                     lambda = -N[[k]]/2+1,
                                     psi = 2*alpha.etak,
                                     chi = sum((Vk-hk)^2/Vk))

        #eta[k]       <- ngme::rGIG( p= -N[[k]]/2+1, a = 2*alpha.etak, b = sum((Vk-hk)^2/Vk))

        start = start + lenghts[k]
        end   = end + lenghts[k+1]
      }
    }
    else{

      post.eta       <- list()
      post.V         <- list()
      SCVI.samples   <- list() #save V matrices for each component
      eta.samples    <- list()

      for(k in 1:ncomp){

        dk     <- d[[k]]
        alpha.etak <- alpha.eta[k]
        hk     <- h[[k]]
        etak   <- eta[k]
        Nk     <- N[[k]]

        if(method == "SVI"){

          ##-------------------------------
          ##     SVI distribution for V's
          ##-------------------------------
          #Read relevant moments from previous iteration first
          if(i>2){
            #Eetam1 <-  fits[[i-1]]$ngvb$post.eta[[k]]$Eetam1
          } else{
            Eetam1 <- 1/etak
          }
          p_V <- -1
          a_V <- Eetam1
          b_V <- dk + hk^2*Eetam1

          ##-------------------------------
          ##     SVI distribution for eta
          ##-------------------------------
          #Get relevant moments from current iteration
          EV   <-  GIGM1(p = p_V, a = a_V, b = b_V)
          EVm1 <-  GIGMm1(p = p_V, a = a_V, b = b_V)
          p_eta <- -length(hk)/2+1
          a_eta <- 2*alpha.etak
          b_eta <- sum(EV-2*hk+hk^2*EVm1)

          ##-----------------------------------------
          ##  Save parameters and relevant moments
          ##-----------------------------------------
          Eeta   <- mGIG(p_eta, a_eta, b_eta, order = 1)   #save gets samples and computes mean by monte carlo
          Eetam1 <- mGIG(p_eta, a_eta, b_eta, order = -1)  #save gets samples and computes mean by monte carlo

          ##----------------------------------------
          ##  Estimate of V to be used in INLA fit
          ##-----------------------------------------
          V[[k]]     <- 1/EVm1
          eta[k]     <- Eeta


          ##----------------------------------------
          ##
          ##----------------------------------------
          post.eta[[k]] <- list(p_eta = p_eta, a_eta = a_eta, b_eta = b_eta, Eeta = Eeta, Eetam1 = Eetam1)
          post.V[[k]]   <- list(p_V = p_V, a_V = a_V, b_V = b_V, EV = EV, EVm1 = EVm1)

        }

        else if(method == "SCVI"){

          ##----------------------------------------
          ##  Estimate of V to be used in INLA fit
          ##-----------------------------------------
          eta.prior <- eta.prior.f(dk,hk,alpha.etak,Nk)
          eta.sample <- sampler.inverseCDF(eta.prior, supp.min = 0, supp.max = 1000, supp.points = 10000, n.samples = n.sampling)

          ntotal <- n.sampling*length(hk)
          p      <- rep(-1,ntotal)
          a      <- rep(1/eta.sample, each=length(hk))
          b      <- c(t(outer(eta.sample, hk, FUN = function(x,y) y^2/x))) + rep(dk,length(eta.sample))
          Vm     <- matrix(1/ngme::rGIG(p,a,b), nrow = n.sampling, ncol = length(hk), byrow = TRUE)
          EVm1   <- colMeans(Vm)

          V[[k]]       <- unname(1/EVm1)
          eta[k]       <- mean(eta.sample)

          SCVI.samples[[k]]  <- 1/Vm
          eta.samples[[k]] <- eta.sample
        }
        else{stop("Method not implemented")}
      }

    }


    ##----------------------------------------
    ##  INLA FIT
    ##----------------------------------------
    fit <- inla.fit.V(V)


    ##----------------------------------------
    ##  Save rest of information
    ##-----------------------------------------

    fit$ngvb <- list(d = d, V = V, eta = eta, time = Sys.time() - start_time)

    if(method == "SVI"){
      fit$ngvb$post.eta <- post.eta
      fit$ngvb$post.V   <- post.V}

    if(method == "SCVI"){  #Save matrix with samples from V
      fit$ngvb$SCVI.samples <- SCVI.samples
      fit$ngvb$eta.samples  <- eta.samples
    }


    fits <- addfit(fits, fit, history)

    #check if Vs change a lot - look at maximum relative variation
    if( (method == "SCVI") || (method == "SVI")){
      old.eta     <- fits@history[[i-1]]$summary.ng$mean
      rate.change <- abs(eta - old.eta)/old.eta
    }

    cat("\n")
    cat(paste("0Oo----------- Iteration ", i ," -----------oO0\n", sep = ""))

    if(verbose){
      cat("Expectation of eta: ")
      cat(round(eta,3))
      cat("\n")

      for(k in 1:ncomp){
        plot(V[[k]]/h[[k]], main = paste("weights for",indexname[k]), ylab = "weights", xlab = "index")
      }
    }


    if( (method == "SCVI") || (method == "SVI") ){
      if(max(rate.change) < stop.rel.change){
        cat(paste("\n 0Oo----------- Convergence achieved -----------oO0\n", sep = ""))
        break}
    }
    if(i == iter){
      cat(paste("\n 0Oo----------- Maximum number of iterations reached -----------oO0\n", sep = ""))
    }

  }


  if(verbose){
    dev.off()
  }

  return(fits)
}

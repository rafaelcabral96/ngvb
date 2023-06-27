
#' @title Diagnostic checks for the latent Gaussianity assumption in R-INLA
#'
#' @description Computes the BF sensitivity measures \eqn{s_0(\mathbf{y})} and  \eqn{I_0(\mathbf{y})} of the non-Gaussianity parameter at
#' the base Gaussian model:
#' \eqn{  \log \pi(\eta|\mathbf{y}) \propto (s_0(\mathbf{y}) - \theta_\eta)\eta  - \frac{I_0(\mathbf{y})}{2} \eta^2 + \mathcal{O}(\eta^2)}, where
#' \eqn{\theta_\eta} is the rate parameter of the exponential PC prior on \eqn{\eta}.
#'
#' @param fit The inla object that fits the LGM.
#' @param Dfunc Function that receives as input the hyperparameter vector \eqn{\boldsymbol{\theta}}
#'     in internal scale and the output is \eqn{\mathbf{D}(\boldsymbol{\theta})}, where \eqn{\mathbf{D}} is the dependency
#'     matrix that specifies the non-Gaussian latent field. If there is more than one latent component
#'     to be extended to non-Gaussianity, this should be a list of functions \eqn{\mathbf{D}_i(\boldsymbol{\theta})},
#'     where \eqn{\mathbf{D}_i(\boldsymbol{\theta})} is the dependency matrix that specifies component i.
#' @param h Predefined constant \eqn{\mathbf{h}} vector that contains the distance between locations
#'     or area of the basis functions. For models defined in discrete space this should be a
#'     vector of ones. If there is more than one latent component to be extended to non-Gaussianity,
#'     this should be a list of vectors \eqn{\mathbf{h_i}}
#'     where \eqn{\mathbf{h_i}} is the predefined constant vector of component i.
#' @param selection List which specifies which model components of the LGM are to be extended
#'     to non-Gaussianity. Same syntax as the argument \code{selection} of
#'     the function \code{inla.posterior.sample}.
#' @param compute.I0    Boolean. If \code{TRUE} compute \eqn{I_0(\mathbf{y}, \hat{\boldsymbol{\gamma}})} (asymptotic variance of the reference
#'  distribution) and p-value.
#' @param compute.predictor Boolean. If \code{TRUE} computes the sensitivity measures for each predictor:
#'      \eqn{\partial_\eta E[z_i|\mathbf{y},\eta]}, where \eqn{z_i} is the i\eqn{^th} element of the linear predictor
#' @param compute.fixed Boolean. If \code{TRUE} computes the sensitivity measures for each fixed effect.
#'      \eqn{\partial_\eta E[z_i|\mathbf{y},\eta]}, where \eqn{z_i} is the i\eqn{^th} fixed effect
#' @param plot Boolean. Generate diagnostic plots.
#' @return A list that for each model component in \code{selection} returns:
#' \itemize{
#'    \item \code{BF.check} The BF sensitivities \eqn{s_0(\mathbf{y})} and \eqn{s_0(\mathbf{y}, \hat{\boldsymbol{\gamma}})}, where
#'    \eqn{\hat{\boldsymbol{\gamma}}} is the posterior mode of the hyperparameters. The variance of the reference distribution
#'    for LGMS with Gaussian response is also shown.
#'    \item \code{BF.index} The contribution of each index to the overall BF sensitivity.
#'    \item \code{sens.fixed.matrix} Sensitivity matrix, containing \eqn{\partial_{\eta_i}E[z_j|\mathbf{y},\eta_i]},
#'    where \eqn{z_j} is the j\eqn{^th} fixed effect and \eqn{\eta_i} is the non-Gaussianity parameter of the i\eqn{^th}
#'    random effect
#'    \item \code{sens.fixed.matrix} Sensitivity matrix, containing \eqn{\partial_{\eta_i}E[z_j|\mathbf{y},\eta_i]},
#'    where \eqn{z_j} is the j\eqn{^th} index of the linear predictor and \eqn{\eta_i} is the non-Gaussianity
#'    parameter of the j\eqn{^th} random effect.
#' }
#'
#' If \code{plot = TRUE}, the BF sensitivities \eqn{d_i(\mathbf{y})} for each index are shown, and also the BF sensitivity at the observed
#' data \eqn{s_0(\mathbf{y}, \hat{\boldsymbol{\gamma}})} is compared to \eqn{s_0(\mathbf{y}^\text{pred}, \hat{\boldsymbol{\gamma}}) \sim N(0,\eqn{I_0(\mathbf{y}), \hat{\boldsymbol{\gamma}}})}. For
#'
#' @examples
#' \dontrun{
#'  #Here we fit an RW1 latent process to the jumpts time series
#'  plot(jumpts)
#'
#'  #Fit LGM with INLA
#'  LGM     <- inla(y ~ -1 + f(x,  model = "rw1"),
#'                  data = jumpts,
#'                  control.compute = list(config = TRUE))
#'
#'  #Fit LnGM with ngvb
#'  check.list <- ng.check(fit = LGM)
#'  }
#' @export
#'
ng.check <- function(fit, Dfunc = NULL, h = NULL, selection = NULL,
                     compute.I0 = FALSE, compute.predictor = FALSE,
                     compute.fixed = TRUE, compute.random = FALSE,
                     plot = TRUE){

  if(is.null(fit$misc$configs)) stop("add control.compute(config = TRUE) argument in INLA call")

  if(is.null(selection)){
    selection <- lapply(fit$summary.random, FUN = function(x) 1: nrow(x))
  }

  ##############################################
  #### Get Dfunc and h #########################
  ##############################################


  comp.names <- names(selection)
  ncomp      <- length(comp.names)

  if(is.null(Dfunc) || is.null(h)){
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
      D.config[[i]] <-  find.D.config(fit, model[i], selection[i], indexname[i], hyperparsel[i,], 0, N[i])
    }
    ##-------------------------------------
    ##    What is h?
    ##-------------------------------------
    h <- list()
    for(i in 1:ncomp){
      h[[i]]          <-  find.h(fit, D.config[[i]])}

    Dfunc <- list()
    for(i in 1:ncomp){
      Dfunc[[i]] <-  D.config[[i]]$Dfunc
    }
  }

  #########################################
  #########################################

  #names of input and output
  random_effect.names <- names(selection)
  fixed_effect.names  <- rownames(fit$summary.fixed)
  if(length(fixed_effect.names) == 0) {compute.fixed <- FALSE}


  sens.fixed.matrix <- NULL
  #get selection list complete
  contents <- fit$misc$configs$contents
  fixed.effect.indeces  <- c()
  predictor.indeces     <- c()

  for(i in 1:length(contents$tag)){
    if(contents$tag[i] == "Predictor"){
      predictor.indeces <- (contents$start[i]):(contents$start[i] + contents$length[i]-1)
    }
    if(contents$tag[i] %in% fixed_effect.names){
      fixed.effect.indeces <- c(fixed.effect.indeces, contents$start[i])
    }
  }

  #sens.fixed.index  <- list()
  if(compute.fixed){ sens.fixed.matrix <- matrix(NA, nrow = length(random_effect.names), ncol = length(fixed_effect.names))}
  #remove Predictor from selection.output
  if(compute.predictor){ sens.pred.matrix <- matrix(NA, nrow = length(random_effect.names), ncol = length(predictor.indeces))}
  if(compute.random){ sens.random <- list()}

  if(compute.fixed | compute.predictor | compute.random){
    for(i in 1:length(random_effect.names)){

      tag.i                 <- which(random_effect.names[i] == contents$tag)
      random.effect.indeces <- (contents$start[tag.i]):(contents$start[tag.i] + contents$length[tag.i] - 1)
      nconfigs <- fit$misc$configs$nconfig

      weights <- numeric(nconfigs)
      for(j in 1:nconfigs){
        weights[j] <- exp(fit$misc$configs$config[[j]]$log.posterior)
      }
      weights <- weights/sum(weights)

      #matrix each row fixed effect, each collum random effect index
      if(compute.fixed){sens.deriv.i <- matrix(0, nrow = length(fixed.effect.indeces), ncol = length(h[[i]]))}
      if(compute.predictor){sens.pred.i  <- matrix(0, nrow = length(predictor.indeces), ncol = length(h[[i]]))}
      if(compute.random){sens.random.i  <- matrix(0, nrow = length(h[[i]]), ncol = length(h[[i]]))}

      for(j in 1:nconfigs){
        mu    <- fit$misc$configs$config[[j]]$mean
        #Sigma <- fit$misc$configs$config[[j]]$Qinv #make sure this works for intrinsic?
        #Sigma <- ultosymmetric(Sigma)              #avoid this - makes it slow... use upper triang version...

        #Sigma1     <- ultosymmetric(fit$misc$configs$config[[j]]$Qinv)
        Sigma     <- solve(ultosymmetric(fit$misc$configs$config[[j]]$Q))

        theta <- fit$misc$configs$config[[j]]$theta
        if(compute.fixed){sens.deriv.i     <- sens.deriv.i +  weights[j]*sens.deriv.i.build(mu, Sigma, Dfunc[[i]](theta), fixed.effect.indeces, random.effect.indeces, h[[i]])}
        if(compute.predictor){sens.pred.i  <- sens.pred.i  +  weights[j]*sens.deriv.i.build(mu, Sigma, Dfunc[[i]](theta), predictor.indeces, random.effect.indeces, h[[i]])}
        if(compute.random){sens.random.i   <- sens.random.i  +  weights[j]*sens.deriv.i.build(mu, Sigma, Dfunc[[i]](theta), random.effect.indeces, random.effect.indeces, h[[i]])}

      }
      if(compute.fixed){
        rownames(sens.deriv.i) <- fixed_effect.names
        sens.fixed.matrix[i,]  <- rowSums(sens.deriv.i)}
      #sens.fixed.index[[i]]  <- sens.deriv.i
      if(compute.predictor){sens.pred.matrix[i,]  <- rowSums(sens.pred.i)}
      if(compute.random){sens.random[[i]]  <- rowSums(sens.random.i)}

    }
    #names(sens.fixed.index) <- random_effect.names
    if(compute.fixed){
      rownames(sens.fixed.matrix) <- random_effect.names
      colnames(sens.fixed.matrix) <- fixed_effect.names}
    if(compute.predictor){rownames(sens.pred.matrix) <- random_effect.names}
  }


  score.list  <- score.mixture.compute(fit, Dfunc, h, selection, compute.I0)


  output <- list(BF.check         = score.list[[1]],
                 BF.index   = score.list[[2]],
                 p.value.index = score.list[[3]],
                 sens.fixed.matrix = sens.fixed.matrix)#, sens.fixed.index = sens.fixed.index)


  if(compute.predictor){
    rownames(sens.pred.matrix) <- random_effect.names
    colnames(sens.pred.matrix) <- 1:length(predictor.indeces)
    output$sens.pred.matrix <- sens.pred.matrix
  }
  if(compute.random){
    names(sens.random) <- random_effect.names
    output$sens.random <- sens.random
  }

  class(output) <- "gaussdiag"
  if(plot) print(plot.gaussdiag(output))


  return(output)
}

#input mu Sigma of the whole latent x, index of the fixed effect and random effect
#output  matrix containing sensititity for each fixed effect and index of the random effect
sens.deriv.i.build <- function(mu, Sigma, D, fixed.effect.indeces, random.effect.indeces, h){

  Szx <- Sigma[fixed.effect.indeces, random.effect.indeces]%*%t(D)
  Sxx <- D%*%Sigma[random.effect.indeces, random.effect.indeces]%*%t(D) #only need diag of this matrix
  muz <- mu[fixed.effect.indeces]
  mux <- (D%*%mu[random.effect.indeces])[,1]

  sens.deriv.i <- matrix(NA, nrow = length(fixed.effect.indeces), ncol = nrow(Sxx))
  for(i in 1:length(fixed.effect.indeces)){

    s12 <- Szx[i,]
    s22 <- diag(Sxx)
    u1  <- muz[i]
    u2  <- mux

    EDx2  <- s22 + u2^2
    EDx4  <- 3*s22^2 + 6*s22*u2^2 + u2^4
    EzDx2 <- 2*s12*u2 + u1*(s22 + u2^2)
    EzDx4 <- 4*s12*(3*s22*u2 + u2^3) + u1*(3*s22^2 + 6*s22*u2^2 + u2^4)

    sens.deriv.i[i,] = (1/(8*h^3))*( 3*h^2*u1 - 6*h*EzDx2 + EzDx4) - u1*(1/(8*h^3))*( 3*h^2 - 6*h*EDx2 + EDx4)


    # for( j in 1:nrow(Sxx)){
    #    s12 <- Szx[i,j]
    #   s22 <- Sxx[j,j]
    #    u1  <- muz[i]
    #   u2  <- mux[j]

    #    EDx2  <- s22 + u2^2
    #   EDx4  <- 3*s22^2 + 6*s22*u2^2 + u2^4

    #    EzDx2 <- 2*s12*u2 + u1*(s22 + u2^2)
    #    EzDx4 <- 4*s12*(3*s22*u2 + u2^3) + u1*(3*s22^2 + 6*s22*u2^2 + u2^4)

    #    sens.deriv.i[i,j] = (1/(8*h[j]^3))*( 3*h[j]^2*u1 - 6*h[j]*EzDx2 + EzDx4) - u1*(1/(8*h[j]^3))*( 3*h[j]^2 - 6*h[j]*EDx2 + EDx4)

    #CzDx2 <- -((-1 + s22)*s22*u1) + 2*s12*u2
    #CzDx4 <- 2*s12*u2 + u1*(s22 + u2^2) - u1*(3*s22^4 + 6*s22^2*u2^2 + u2^4)
    #sens.deriv.i[i,j] = (1/(8*h[j]^3))*(- 6*h[j]*CzDx2 + CzDx4)
    # }
  }
  return(sens.deriv.i)
}


score.mixture.compute <- function(fit, Dfunc, h, selection, compute.I0){

  n.config <- fit$misc$configs$nconfig          #number of hyperparameter configurations

  mean.list  <- list()
  Qinv.list  <- list()
  theta.list <- list()
  delta.vec  <- c()

  for(i in 1:n.config){
    mean.list[[i]]  <- fit$misc$configs$config[[i]]$mean
    Qinv.list[[i]]  <- fit$misc$configs$config[[i]]$Qinv #avoid this use Qinv try solve here...
    if(n.config == 1){
      if(fit$.args$control.inla$int.strategy=="eb"){
        theta.list[[i]]  <- fit$mode$theta
      }
      else{
        theta.list[[i]] <- fit$.args$control.mode$theta}# fit$mode$theta}
    }
    else{theta.list[[i]] <- fit$misc$configs$config[[i]]$theta}
    delta.vec[i]    <- exp(fit$misc$configs$config[[i]]$log.posterior)
  }
  delta.vec  <- delta.vec/sum(delta.vec)

  comp.names <- names(selection)
  n.comp     <- length(comp.names)



  score.list          <- replicate(n.comp,list())         #score S_0(y,gamma)
  I0.list             <- replicate(n.comp,list())         #information I_0(y,gamma)
  variance.ref.list   <- replicate(n.comp,list())         #variance of the reference distribution for each gamma
  p.value.global.list <- replicate(n.comp,list())         #pvalue(y,gamma)
  p.value.index.list  <- replicate(n.comp,list())
  p.value.index       <- list()
  score.index         <- list()
  score               <- numeric(n.comp)
  s0.mode             <- numeric(n.comp)
  var.ref.mode        <- numeric(n.comp)
  I0                  <- numeric(n.comp)
  p.value.global      <- numeric(n.comp)


  if(compute.I0){
    I0      <- numeric(n.comp)   #I0(y)
    I0.list <- list()            #I0(y,gamma)
  }

  for(j in 1:n.comp){
    hj      <- h[[j]]

    sel.tag <- fit$misc$configs$contents$tag == comp.names[j]
    init    <- fit$misc$configs$contents$start[sel.tag]
    end     <- init + fit$misc$configs$contents$length[sel.tag] -1
    sel     <- init:end

    #computing mean and variances od [DX]---------------------------
    for(k in 1:n.config){
      D <- Dfunc[[j]](theta.list[[k]])

      mu.Dx.k    <- (D%*%mean.list[[k]][sel])[,1]

      Sigma.Dx.k <- (D%*%ultosymmetric(Qinv.list[[k]][sel,sel])%*%t(D)) #avoid this - makes it slow... use upper triang version...
      var.Dx.k   <- diag(Sigma.Dx.k)

      d.2.k <- ( mu.Dx.k^2 + var.Dx.k )
      d.4.k <- ( mu.Dx.k^4 + 6*mu.Dx.k^2*var.Dx.k + 3*var.Dx.k^2 )
      score.list[[j]][[k]] <- ( 3*hj^2 - 6*hj*d.2.k + d.4.k ) / ( 8*hj^3 )

      if(compute.I0){
        d.6 <- mu.Dx.k^6 + 15*mu.Dx.k^4*var.Dx.k + 45*mu.Dx.k^2*var.Dx.k^2 + 15*var.Dx.k^3

        weird1 <- weird1compute(mu.Dx.k,var.Dx.k,hj)
        weird2 <- weird2compute(mu.Dx.k,as.matrix(Sigma.Dx.k),hj)
        weird  <- sum(weird1) +  2*weird2
        Varf   <- weird - sum(score.list[[j]][[k]])^2

        I0.list[[j]][[k]]             <- sum( ( 3*hj^3 + 3*hj^2*d.2.k -6*hj*d.4.k + d.6.k ) / ( 8*hj^5 ) )  - Varf
        p.value.global.list[[j]][[k]] <-  1 - pnorm(sum(score.list[[j]][[k]]), sd = sqrt(I0.list[[j]][[k]]))

      }

      #compute variance of the reference distribution
      if((!compute.I0) && fit$.args$family == "gaussian"){
        Gamma                          <- -Sigma.Dx.k + diag(h[[j]])
        p.value.index.list[[j]][[k]]   <-  pvalue.d(score.list[[j]][[k]], h[[j]], diag(Gamma))

        variance.ref.list[[j]][[k]]   <-  (3/8)*(t(h[[j]]^(-3))%*%(Gamma^4)%*%(h[[j]]^(-3)))[1,1]
        p.value.global.list[[j]][[k]] <-  1 - pnorm(sum(score.list[[j]][[k]]), sd = sqrt(variance.ref.list[[j]][[k]]))

       }

    }

    #Now do weighted average - score
    score.index[[j]]    <- custom.weighted.mean(score.list[[j]], delta.vec)
    score[j]            <- sum(score.index[[j]])

    #Now do weighted average - p.value global
    if(compute.I0 || fit$.args$family == "gaussian"){
      p.value.global[j] <- custom.weighted.mean(p.value.global.list[[j]], delta.vec)
    }

    #Now do weighted average - p.value single loc
    if(fit$.args$family == "gaussian"){
      p.value.index[[j]]   <- custom.weighted.mean(p.value.index.list[[j]], delta.vec)
    }

    #s0 mode
    s0.mode[j] = lapply(score.list[[j]], sum)[[which.max(delta.vec)]]
    if(fit$.args$family == "gaussian"){
      var.ref.mode[j] = variance.ref.list[[j]][[which.max(delta.vec)]]
    }
    else if (compute.I0){
      var.ref.mode[j] = I0.list[[j]][[which.max(delta.vec)]]
    }

  }

  names(score.index)    <- comp.names
  names(p.value.index)  <- comp.names
  names(score)          <- comp.names
  score.dataframe       <- data.frame(s0 = score, s0.mode = s0.mode, var.ref.mode = var.ref.mode)

  if(compute.I0 || fit$.args$family == "gaussian"){
    score.dataframe$p.value      = p.value.global
  }

  #names(p.value.index)     <- comp.names
  #names(score.list)        <- comp.names
  #names(variance.ref.list) <- comp.names
  #names(I0.list)           <- comp.names

  #misc = list(score.list = score.list, variance.ref.list = variance.ref.list, I0.list = I0.list, delta.vec = delta.vec)

  return(list(BF.check = score.dataframe, BF.index = score.index, p.value.index))
}

ultosymmetric=function(m){
  m = m + t(m) - diag(diag(m))
  return (m)}

cov.2m <- function(input, output){
  n <- ncol(input)
  m <- ncol(output)

  x <-foreach(i=1:n, .combine='cbind') %:%
    foreach(j=1:m, .combine='c') %do% {
      #cov(input[,i], output[,j])
      mean(input[,i]*(output[,j]-mean(output[,j])))}

  return(x)
}

pvalue.d <- function(x, h, g){ #paralalelize in c++ make more efficient
  y <- numeric(length(x))
  for(i in 1:length(x)){
    if(x[i]>(3*g[i]^2)/(8*h[i]^3)){
      y[i] = 2*(1 - pnorm(sqrt(3+sqrt(6+8*h[i]^3*g[i]^(-2)*x[i]))))
    }
    else if(x[i] > -(3*g[i]^2)/(4*h[i]^3)){
      y[i] = 1 +2*( pnorm(sqrt(3-sqrt(6+8*h[i]^3*g[i]^(-2)*x[i]))) - pnorm(sqrt(3+sqrt(6+8*h[i]^3*g[i]^(-2)*x[i]))))
    }
    else{
      y[i] = 1
    }
  }
  return(y)
}


custom.weighted.mean <- function(list, w){
  n <- length(w)
  if(n == 1){return(list[[1]])}
  else{
    out <- w[1]*list[[1]]
    for(j in 1:n){
      out <- out + w[j]*list[[j]]
    }
  }
  return(out)
}

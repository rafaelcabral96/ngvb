if(0){


  #' @title Diagnostic checks for the latent Gaussianity assumption in R-INLA
  #'
  #' @description Computes the observed Fisher's score \eqn{s} and information \eqn{I} of the non-Gaussianity parameter at
  #' the base Gaussian model:
  #' \eqn{  \log \pi(\eta|\mathbf{y}) \propto (s - \theta_\eta)\eta  - \frac{I}{2} \eta^2 + \mathcal{O}(\eta^2)}, where
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
  #' @param compute.fisher    Boolean. If \code{TRUE} compute the observed fisher information and p-value.
  #' @param compute.predictor Boolean. If \code{TRUE} computes the sensitivity measures for each predictor:
  #'      \eqn{\partial_\eta E[z_i|\mathbf{y},\eta]}, where \eqn{z_i} is the i\eqn{^th} element of the linear predictor
  #' @param compute.fixed Boolean. If \code{TRUE} computes the sensitivity measures for each fixed effect.
  #'      \eqn{\partial_\eta E[z_i|\mathbf{y},\eta]}, where \eqn{z_i} is the i\eqn{^th} fixed effect
  #' @param method Either "mixture" or "sampling". Compute observed Fisher's score, information and sensitivity measures
  #'      analytically using the Gaussian mixture approximation or from the posterior samples obtained by
  #'      \code{inla.posterior.sample}.
  #' @param n.sample Integer. Number of samples to use in \code{inla.posterior.sample}.
  #' @param plot Boolean. Generate diagnostic plots.
  #' @return A list that for each model component in \code{selection} returns:
  #' \itemize{
  #'    \item \code{score} The observed Fisher's score \eqn{s} and information \eqn{I}, along with the Bayesian p-value
  #'    obtained from the asymptotic reference distribution \eqn{s^{\text{ref}} \sim N(0,I)}.
  #'    \item \code{score.index} The contribution of each index to the overall score.
  #'    \item \code{sens.fixed.matrix} Sensitivity matrix, containing \eqn{\partial_{\eta_i}E[z_j|\mathbf{y},\eta_i]},
  #'    where \eqn{z_j} is the j\eqn{^th} fixed effect and \eqn{\eta_i} is the non-Gaussianity parameter of the i\eqn{^th}
  #'    random effect
  #'    \item \code{sens.fixed.matrix} Sensitivity matrix, containing \eqn{\partial_{\eta_i}E[z_j|\mathbf{y},\eta_i]},
  #'    where \eqn{z_j} is the j\eqn{^th} index of the linear predictor and \eqn{\eta_i} is the non-Gaussianity
  #'    parameter of the j\eqn{^th} random effect.
  #' }
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
  #'  score.list <- ng.score(fit = LGM)
  #'  }
  #' @export
  #'
  ng.score <- function(fit, Dfunc = NULL, h = NULL, selection = NULL,
                       compute.fisher = FALSE, compute.predictor = FALSE, compute.fixed = TRUE, compute.random = FALSE,
                       method = "mixture", n.sample = 5000, plot = TRUE){

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

    if(method == "sampling"){

      sens.fixed.matrix <- NULL

      #get selection list complete
      contents <- fit$misc$configs$contents
      selection.all <- list()
      k <- 1
      for(i in 1:length(contents$tag)){
        if(compute.predictor & contents$tag[i] == "Predictor"){
          selection.all[[k]] <- 1:contents$length[i]
          names(selection.all)[k] <- contents$tag[i]
          k <- k + 1
        }
        if(contents$tag[i] %in% c(random_effect.names, fixed_effect.names)){
          selection.all[[k]] <- 1:contents$length[i]
          names(selection.all)[k] <- contents$tag[i]
          k <- k + 1
        }
      }

      samples = inla.posterior.sample(n = n.sample, fit, selection = selection.all,
                                      intern = TRUE, use.improved.mean = TRUE, skew.corr = TRUE,
                                      add.names = TRUE, seed = 0L, num.threads = 0,
                                      verbose = FALSE)

      ls <- c(0,cumsum(lengths(selection.all)))

      out <- foreach(j = 1:n.sample, .combine = 'rbind', .packages = "Matrix") %do% {
        hyperpar  <- samples[[j]]$hyperpar
        sample <- samples[[j]]

        input.list <- list()
        fisher.list <- list()
        random.list <- list()
        for(k in 1:length(random_effect.names)){
          index.k         <- match(random_effect.names[k], names(ls))
          sample.k        <- sample$latent[(ls[index.k-1]+1):ls[index.k]]
          Lambda.k        <- (Dfunc[[k]](hyperpar)%*%sample.k)[,1]
          input.list[[k]] <- 1/(8*h[[k]]^3)*(3*h[[k]]^2 - 6*h[[k]]*Lambda.k^2 + Lambda.k^4)

          if(compute.random){ random.list[[k]] <- sample.k}
          if(compute.fisher) {fisher.list[[k]] <- 1/(8*h[[k]]^5)*(3*h[[k]]^3 + 3*h[[k]]^2*Lambda.k^2 - 6*h[[k]]*Lambda.k^4 + Lambda.k^6)}
        }
        names(input.list) <- random_effect.names
        if(compute.fisher) {names(fisher.list) <- random_effect.names}
        if(compute.random) {names(random.list) <- random_effect.names}

        output.list <- list()
        k <- 1

        if(compute.fixed){
          for(k in 1:length(fixed_effect.names)){
            index.k         <- match(fixed_effect.names[k], names(ls))
            output.list[[k]]  <- sample$latent[(ls[index.k-1]+1):ls[index.k]]
          }
          names(output.list) <- fixed_effect.names
          k <- k + 1
        }

        if(compute.predictor){
          index.k                  <- match("Predictor", names(ls))
          output.list[[k]]         <- sample$latent[(ls[index.k-1]+1):ls[index.k]]
          names(output.list)[[k]]  <-  "Predictor"
        }


        return(list(input = input.list, output = output.list, fisher = fisher.list, random = random.list))}


      input.list  <- out[1:n.sample]
      output.list <- out[(n.sample+1):(2*n.sample)]
      fisher.list <- out[(2*n.sample+1):(3*n.sample)]
      random.list <- out[(3*n.sample+1):(4*n.sample)]

      score <- c()
      score.index <- list()
      fisher <- c()


      #sens.fixed.index <- list()
      if(compute.fixed){sens.fixed.matrix <- matrix(NA, nrow = length(random_effect.names), ncol = length(fixed_effect.names))}

      #remove Predictor from selection.output
      if(compute.predictor){
        sens.pred.matrix <- matrix(NA, nrow = length(random_effect.names), ncol = length(selection.all$Predictor))
        pred.matrix  <- do.call(rbind, lapply(output.list, FUN = "[[", "Predictor"))
      }
      if(compute.random){sens.random  <- list()}

      for(i in 1:length(random_effect.names)){
        input.matrix  <- do.call(rbind, lapply(input.list, FUN = "[[", random_effect.names[i]))
        score.index[[i]]   <- colMeans(input.matrix)
        score <- c(score, sum(score.index[[i]]))

        if(compute.fixed){
          cov.i <- matrix(NA, nrow = length(fixed_effect.names), ncol = ncol(input.matrix))
          for(j in 1:length(fixed_effect.names)){
            output.vector <- do.call(c,  lapply(output.list,   FUN = "[[", fixed_effect.names[j]))
            cov.ij <- cov.2m(input.matrix, as.matrix(output.vector))
            sens.fixed.matrix[i,j] <- sum(cov.ij)
            cov.i[j,] <- cov.ij
          }
          rownames(cov.i) <- fixed_effect.names
          #sens.fixed.index[[i]] <- cov.i
        }

        if(compute.fisher){
          fisher.matrix <- do.call(rbind, lapply(fisher.list, FUN = "[[", random_effect.names[i]))
          fisher <- c(fisher, sum( colMeans(fisher.matrix - input.matrix^2)) +  score[i]^2)
        }
        if(compute.predictor){sens.pred.matrix[i,] <- rowSums(cov.2m(input.matrix, pred.matrix))}
        if(compute.random){
          random.matrix <-  do.call(rbind, lapply(random.list, FUN = "[[", random_effect.names[i]))
          sens.random[[i]]    <- rowSums(cov.2m(input.matrix, random.matrix))}
      }

      names(score)                <- random_effect.names
      names(score.index)          <- random_effect.names
      if(compute.fixed){
        #names(sens.fixed.index)     <- random_effect.names
        rownames(sens.fixed.matrix) <- random_effect.names
        colnames(sens.fixed.matrix) <- fixed_effect.names
      }

      if(compute.fisher) {
        names(fisher) <- random_effect.names
        score <- data.frame( score = score, fisher = fisher, p.value = 1 - pnorm(score, sd = sqrt(fisher)))
      }else{ score <- data.frame( score = score ) }
      output <- list(score = score, score.index = score.index, sens.fixed.matrix = sens.fixed.matrix)#, sens.fixed.index = sens.fixed.index)

    }
    else if(method == "mixture"){

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


      score.list  <- score.mixture.compute(fit, Dfunc, h, selection, compute.fisher)


      output <- list(score = score.list[[1]],
                     score.index = score.list[[2]],
                     p.value.index = score.list[[3]],
                     dist.index = score.list[[4]],
                     sens.fixed.matrix = sens.fixed.matrix)#, sens.fixed.index = sens.fixed.index)

    }

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


  score.mixture.compute <- function(fit, Dfunc, h, selection, compute.fisher){

    n.config <- fit$misc$configs$nconfig          #number of hyperparameter configurations

    mean.list  <- list()
    Qinv.list  <- list()
    theta.list <- list()
    delta.vec  <- c()

    for(i in 1:n.config){
      mean.list[[i]]  <- fit$misc$configs$config[[i]]$mean
      Qinv.list[[i]]  <- fit$misc$configs$config[[i]]$Qinv
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

    score        <- numeric(n.comp)
    variance.ref <- rep(0,n.comp) #variance of the reference distribution
    score.index <- list()
    p.value.index <- list()
    dist.index  <- list()

    if(compute.fisher){
      fisher      <- numeric(n.comp)
      Varf        <- numeric(n.comp) #KLD(post flex || post base ) = Varf eta^2
    }

    for(j in 1:n.comp){
      hj      <- h[[j]]

      sel.tag <- fit$misc$configs$contents$tag == comp.names[j]
      init    <- fit$misc$configs$contents$start[sel.tag]
      end     <- init + fit$misc$configs$contents$length[sel.tag] -1
      sel     <- init:end


      #computing mean and variances od [DX]---------------------------
      d.2 <- 0
      d.4 <- 0
      d.6 <- 0
      weird <- 0
      KLDj  <- 0
      p.value.index.j <- 0

      for(k in 1:n.config){
        D <- Dfunc[[j]](theta.list[[k]])

        mu.Dx.k    <- (D%*%mean.list[[k]][sel])[,1]

        Sigma.Dx.k <- (D%*%ultosymmetric(Qinv.list[[k]][sel,sel])%*%t(D)) #avoid this - makes it slow... use upper triang version...
        var.Dx.k   <- diag(Sigma.Dx.k)

        d.2.k <- ( mu.Dx.k^2 + var.Dx.k )
        d.4.k <- ( mu.Dx.k^4 + 6*mu.Dx.k^2*var.Dx.k + 3*var.Dx.k^2 )

        d.2 <- d.2 + delta.vec[k]*d.2.k
        d.4 <- d.4 + delta.vec[k]*d.4.k

        if(compute.fisher){
          d.6 <- d.6 + delta.vec[k]*( mu.Dx.k^6 + 15*mu.Dx.k^4*var.Dx.k + 45*mu.Dx.k^2*var.Dx.k^2 + 15*var.Dx.k^3 )

          weird1 <- weird1compute(mu.Dx.k,var.Dx.k,hj)
          weird2 <- weird2compute(mu.Dx.k,as.matrix(Sigma.Dx.k),hj)
          weird  <- weird + delta.vec[k]*( sum(weird1) +  2*weird2 )
          KLDj   <- KLDj + delta.vec[k]*weird1 #increase in KLD per point
        }

        #compute variance of the reference distribution
        if(fit$.args$family == "gaussian"){
          Gamma           <- -Sigma.Dx.k + diag(h[[j]])
          variance.ref[j] <-  variance.ref[j] + delta.vec[k]*3/8*(t(h[[j]]^(-3))%*%(Gamma^4)%*%(h[[j]]^(-3)))[1,1]
          s0.k            <- ( 3*hj^2 - 6*hj*d.2.k + d.4.k ) / ( 8*hj^3 )  ##make this code more efficient, I am repeating same computation several times
          p.value.index.j <-  p.value.index.j + delta.vec[k]*pvalue.d(s0.k, h[[j]], diag(Gamma))
        }

      }


      deriv              <- ( 3*hj^2 - 6*hj*d.2 + d.4 ) / ( 8*hj^3 )
      score[j]           <- sum( deriv )
      score.index[[j]]   <- deriv

      if(fit$.args$family == "gaussian"){
        p.value.index[[j]] <- p.value.index.j
      }

      if(compute.fisher){
        Varf[j]        <- weird - score[j]^2
        fisher[j]      <- sum( ( 3*hj^3 + 3*hj^2*d.2 -6*hj*d.4 + d.6 ) / ( 8*hj^5 ) )  - Varf[j]
        dist.index[[j]]<- sqrt(2*KLDj)/sqrt((3/8)*hj^(-2)) #increase of distance in posterior / increase of distance in prior for each point...
      }

    }

    names(score.index)    <- comp.names
    names(score)          <- comp.names

    score.dataframe = data.frame(score = score)

    if(compute.fisher){
      names(dist.index) <- comp.names

      dist.post <- sqrt(2*Varf)  #distance increase in posterior
      dist.pert <- sapply(h, function(x) sqrt((3/8)*sum(x^(-2)))) #distance increase in prior

      score.dataframe$fisher  = fisher
      score.dataframe$dist    = dist.post/dist.pert
      score.dataframe$p.value = p.value =  1 - pnorm(score, sd = sqrt(fisher))
    }

    if(fit$.args$family == "gaussian"){
      names(variance.ref)   <- comp.names
      names(p.value.index)  <- comp.names

      score.dataframe$variance.ref = variance.ref
      score.dataframe$p.value = p.value =  1 - pnorm(score, sd = sqrt(variance.ref))
    }

    return(list(score =  score.dataframe, score.index, p.value.index, dist.index))
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

  eoutliers <- function(y, fit, n.sample){

    if(fit$.args$family != "gaussian") stop("family must be Gaussian")
    if(is.null(fit$misc$configs)) stop("add control.compute(config = TRUE) argument in INLA call")

    N <- fit$size.linear.predictor$N

    samples  <- inla.posterior.sample(n.sample, fit, selection = list("Predictor" = 1:N))

    s <- rep(0,length(y))
    for(i in 1:n.sample){
      lp <- samples[[i]]$latent[,1]
      sigma_y <- 1/sqrt(samples[[i]]$hyperpar[1])
      r       <- (y - lp)/sigma_y
      s       <- s + (3 -6*r^2+r^4)/8
    }
    s <- s/n.sample

    return(list(score = sum(s), score.index = s))
  }


  #' @title Diagnostic checks for the latent Gaussianity assumption in R-INLA
  #'
  #' @description Computes the observed Fisher's score \eqn{s} and information \eqn{I} of the non-Gaussianity parameter at
  #' the base Gaussian model:
  #' \eqn{  \log \pi(\eta|\mathbf{y}) \propto (s - \theta_\eta)\eta  - \frac{I}{2} \eta^2 + \mathcal{O}(\eta^2)}, where
  #' \eqn{\theta_\eta} is the rate parameter of the exponential PC prior on \eqn{\eta}.
  #'
  #' @param fit The inla object that fits the LGM.
  #' @param manual.configs Not necessary if \code{fit} is provided. List containing:
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
  #' @param ranking.only Boolean. If \code{TRUE} compute only the score \eqn{s}.
  #' @param plot Boolean. Generate diagnostic plots.
  #' @return A list that for each model component in \code{selection} returns:
  #' \itemize{
  #'    \item \code{summary} The observed Fisher's score \eqn{s} and information \eqn{I}, along with the scale of the
  #'    asymptotic reference distribution \eqn{s^{\text{ref}} \sim N(0,I)} and the Bayesian p-value.
  #'    \item \code{deriv} The contribution of each index to the overall score.
  #'  }
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
  #'  score.list <- ng.score(fit = LGM, selection = list(x=1:100))
  #'
  #'  #Visual diagnostic
  #'  plot.score(score.list, largest = 2)[[3]]
  #'  }
  #' @export
  #'
  ng.score2 <- function(fit, Dfunc = NULL, h = NULL, selection = NULL, ranking.only = FALSE){

    if(is.null(fit$misc$configs)) stop("add control.compute(config = TRUE) argument in INLA call")

    if(is.null(selection)){
      selection <- lapply(fit$summary.random, FUN = function(x) 1: nrow(x))
    }

    ##-------------------------------------
    ##    HOW MANY COMPONENTS TO CHECK
    ##-------------------------------------

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

    diagnostic <- list()

    I0 <- NA

    #extracting relevant configurations-------------------------------
    n.config <- fit$misc$configs$nconfig          #number of hyperparameter configurations

    mean.list  <- list()
    Qinv.list  <- list()
    theta.list <- list()
    delta.vec  <- c()

    for(i in 1:n.config){
      mean.list[[i]]  <- fit$misc$configs$config[[i]]$mean
      Qinv.list[[i]]  <- fit$misc$configs$config[[i]]$Qinv
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


    for(j in 1:n.comp){
      hj      <- h[[j]]

      sel.tag <- fit$misc$configs$contents$tag == comp.names[j]
      init    <- fit$misc$configs$contents$start[sel.tag]
      end     <- init + fit$misc$configs$contents$length[sel.tag] -1
      sel     <- init:end


      #computing mean and variances od [DX]---------------------------
      d.2 <- 0
      d.4 <- 0
      d.6 <- 0
      weird <- 0

      for(k in 1:n.config){
        D <- Dfunc[[j]](theta.list[[k]])

        mu.Dx.k    <- (D%*%mean.list[[k]][sel])[,1]

        Sigma.Dx.k <- (D%*%ultosymmetric(Qinv.list[[k]][sel,sel])%*%t(D)) #avoid this - makes it slow... use upper triang version...
        var.Dx.k   <- diag(Sigma.Dx.k)

        d.2 <- d.2 + delta.vec[k]*( mu.Dx.k^2 + var.Dx.k )
        d.4 <- d.4 + delta.vec[k]*( mu.Dx.k^4 + 6*mu.Dx.k^2*var.Dx.k + 3*var.Dx.k^2 )

        if(!ranking.only){
          d.6 <- d.6 + delta.vec[k]*( mu.Dx.k^6 + 15*mu.Dx.k^4*var.Dx.k + 45*mu.Dx.k^2*var.Dx.k^2 + 15*var.Dx.k^3 )

          weird1 <- weird1compute(mu.Dx.k,var.Dx.k,hj)
          weird2 <- weird2compute(mu.Dx.k,as.matrix(Sigma.Dx.k),hj)
          weird  <- weird + delta.vec[k]*( weird1 +  2*weird2 )
        }

      }


      deriv <- ( 3*hj^2 - 6*hj*d.2 + d.4 ) / ( 8*hj^3 )
      s0    <- sum( deriv )

      if(!ranking.only){
        I0    <- sum( ( 3*hj^3 + 3*hj^2*d.2 -6*hj*d.4 + d.6 ) / ( 8*hj^5 ) )  + s0^2 - weird
      }

      scale          <- sqrt(I0)
      pValue         <- 1-pnorm(s0, mean = 0, sd = scale)
      summary        <- c( s0, I0, scale, pValue)
      names(summary) <- c("Score", "Observed Fisher Info.", "scale", "p-value")

      #convert numeric intro matrix to have same output type

      diagnostic[[j]] <- list(summary = summary, deriv = deriv)




    }

    names(diagnostic) <- comp.names

    #class(diagnostic) <- "scorelist"


    #if(plot) print(plot(diagnostic))

    return(diagnostic)

  }

  }

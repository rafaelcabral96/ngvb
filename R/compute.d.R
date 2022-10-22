compute.d.sampling <- function(fit, D.config){

  ncomp       <- length(D.config)
  d           <- list()
  for(i in 1:ncomp){
    n.sampling  <- D.config[[i]]$n.sampling
    selection   <- D.config[[i]]$selection
    Dfunc       <- D.config[[i]]$Dfunc

    samples <- inla.posterior.sample(n = n.sampling, fit, selection = selection,
                                     intern = TRUE, use.improved.mean = TRUE, skew.corr = TRUE,
                                     add.names = TRUE, seed = 0L, num.threads = 0,
                                     verbose = FALSE)


    Lambda <- foreach(j = 1:n.sampling, .combine = 'rbind', .packages = "Matrix") %do% {
      hyperpar  <- samples[[j]]$hyperpar
      postx     <- c(samples[[j]]$latent)
      (Dfunc(hyperpar)%*%postx)[,1]}

    d[[i]] <- colMeans(Lambda^2)
  }

  return(d)
}

compute.d.configs <- function(fit, D.config){


 # a <- matrix(rnorm(3^2), nrow=3)
#  print(a)
#  print(t(a))

#  print(showMethods(t))
#  a <- AR1.matrix(100,0.9)
#  print(a)
#  print(dgCMatrix::t(a))


  ncomp       <- length(D.config)
  d           <- list()
  for(i in 1:ncomp){
    indexname   <- D.config[[i]]$indexname
    Dfunc       <- D.config[[i]]$Dfunc
    h           <- D.config[[i]]$h

    nconfigs <- length(fit$misc$configs$config)
    ll       <- numeric(nconfigs)                 #ll for each weight
    dconfigs <- matrix(NA, nrow = nconfigs, ncol = length(h))

    #matrix Q has high dimension predictor + x + coefficients, etc. you need to pick relevant position
    field.pos = which(fit$misc$configs$contents$tag == indexname)
    matrix.divisions = c(0,cumsum(fit$misc$configs$contents$length[-1]))
    matrix.start = matrix.divisions[field.pos-1]+1
    matrix.end   = matrix.divisions[field.pos]

    dmatrix   <- matrix(NA, nrow = nconfigs, ncol = length(h))
    for(j in 1:nconfigs){
      configi <- fit$misc$configs$config[[j]]
      m       <- configi$improved.mean[matrix.start:matrix.end]
      #Sigma   <- chol2inv(configi$Q[matrix.start:matrix.end, matrix.start:matrix.end]) #this is not the precision matrix of the conditional posterior
      Sigma    <- configi$Qinv[matrix.start:matrix.end, matrix.start:matrix.end] #not Sigma this is diagonal
      Di       <- Dfunc(configi$theta)     #Attention configi$theta is in internal scale
      dmatrix[j,] <- ((Di%*%m)^2)[,1] +  diag(Di%*%Sigma%*%t(Di))
      ll[j]   <- configi$log.posterior
    }
    weights = exp(ll)/sum(exp(ll))
    di <- rep(0, length(h))
    for(j in 1:nconfigs){
      di <- di + dmatrix[j,]*weights[j]
    }
    d[[i]] <- di
  }

  return(d)
}


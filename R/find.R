find.h <- function(fit, D.config){
  model        <- D.config$model
  N            <- D.config$N
  indexname    <- D.config$indexname
  Dfunc        <- D.config$Dfunc

  if( model == "rw1" ){
    x  <- fit$.args$data[[indexname]]
    h  <- diff(x)
  } else if ( model == "rw2"){
    h  <- rep(1, N-2)
  } else if( any(model == c("ar1","iid")) ){
    h  <- rep(1, N)
  }
  return(h)
}

find.D.config <- function(fit, model, selection, indexname, hyperparsel, n.sampling, N){

  Dfunc             <- Dfuncsel(model, N, D.config = list(hyperparsel = hyperparsel))
  D                 <- Dfunc(rep(0,length(hyperparsel)))

  D.config          <- list(hyperparsel = hyperparsel, D = D, Dfunc = Dfunc,
                            selection = selection, indexname = indexname,
                            n.sampling = n.sampling, N = N, model = model)

  if(model == "ar1"){
    D.config$pcprior <- ar1.get.prior(fit, indexname)
  }

  return(D.config)

}

#return inla.fit.V(V, theta) object
#If theta = NULL, then theta is a model parameter
#If theta is assigned then it is a fixed parameter
find.inla.fit.V <- function(fit, N, model, term_string, ncomp, D.config){

  ##----------------------------------------
  ##  New INLA fit arguments
  ##----------------------------------------
  arglist <- fit$.args
  str.formula <- gsub(" ", "", paste(format(fit$.args$formula), collapse = ""), fixed = TRUE)
  for(i in 1:ncomp){
    modeli   <- model[[i]]
    termiold <- gsub(" ", "", term_string[i], fixed = TRUE)

    if(modeli == "ar1") { #use ar1 rgeneric
      terminew <- gsub(paste("\"",model[i],"\"", sep = ""), paste("rgeneric",i,sep=""),termiold) #model = "rw1"
      terminew <- sub(",hyper.*)", ")", terminew)   #remove hyper = prior.list from formula so you can use cgeneric
    }
    else if(modeli == "iid") { #use iid scales
      terminew <- gsub(")", paste0(", scale = 1/V[[",i,"]] )", sep=""), termiold) #model = "rw1"
    }
    else{  #use generic0
      terminew <- gsub(paste("\"",model[i],"\"", sep = ""), paste("\"generic0\", Cmatrix = VINIG.Q.",i,", rankdef = VINIG.rank.def.",i, sep=""),termiold) #model = "rw1"
    }
    str.formula <- gsub(termiold, terminew, str.formula, fixed = T)
  }
  arglist$formula <- formula(str.formula)

  ##---------------------------------------
  ## create rgneric model
  ##---------------------------------------

  #Use generic0 or cgeneric for componenti?
  #This chunk of code can be moved to D.config.generate part
  use.generic0 = rep(FALSE, ncomp)

  for(i in 1:ncomp){
    if(model[i] == "rw1" || model[i] == "rw2"){
      #D.config[[i]]$log.det.D       <- det.mp(D.config[[i]]$D) #IMPROVE THIS! now it it now needed
      if(model[i] == "rw1")   D.config[[i]]$rankdef = 1
      if(model[i] == "rw2")   D.config[[i]]$rankdef = 2
      assign(paste('VINIG.rank.def.',i,sep=""), D.config[[i]]$rankdef, envir = .GlobalEnv)
      use.generic0[i] = TRUE
      print("Warning: You need to add log(|D|) to mll: check inla.doc(generic0)")
    }
  }


  inla.fit.V <- function(V, theta = NULL){

    for(i in 1:ncomp){

      #USE generic0 or cgeneric?
      if(use.generic0[i]){ #If so you only need to define matrices Q and rankdefs
        D <- D.config[[i]]$D

        assign(paste('VINIG.Q.',i,sep=""), as(t(D)%*%Diagonal(length(V[[i]]),1/V[[i]])%*%D, "dgTMatrix"), envir = .GlobalEnv)

      }else if( model[i] == "iid") {
        assign("V", V, envir = .GlobalEnv)

      }else if (model[i] == "ar1"){ #If so you only need to define matrices Q and rankdefs


        pcprior <- D.config[[i]]$pcprior

        folder  <- paste0(system.file(package = 'ngvb'),'/cfiles/')

        cmodel  <- inla.cgeneric.define(model = "inla_cgeneric_ar1_model",
                                        shlib = paste0(folder,"cgeneric-ngvb.so"),
                                        n     = N[[i]],
                                        V     = V[[i]],
                                        U1    = pcprior$U1,
                                        alpha1= pcprior$alpha1,
                                        U2    = pcprior$U2,
                                        alpha2= pcprior$alpha2)

        assign(paste('rgeneric',i,sep=""), cmodel, envir = .GlobalEnv)

      }

      #needed for sampling from x for a fixed value of theta
      if(!is.null(theta)){
        arglist$control.mode = list(theta = theta, fixed = TRUE)
      }

      fit       <- do.call(inla, arglist)
      return(fit)
    }

    return(inla.fit.V)
  }
}

Dfuncsel <- function(model, N, D.config = NULL){
  hyperparsel = D.config$hyperparsel
  if(model == "rw1"){ #model with other boundary conditions not implemented
    D     <- RW1.matrix(N)
    Dfunc <- function(theta){
      prec <- exp(theta[hyperparsel])
      return(sqrt(prec)*D)
    }
  }

  else if(model == "rw2"){ # this is it in discrete case, for continuous case look at H paper
    D     <- RW2.matrix(N)
    Dfunc <- function(theta){
      prec <- exp(theta[hyperparsel])
      return(sqrt(prec)*D)
    }
  }

  else if(model == "ar1"){
    Dfunc <- function(theta){
      theta           <- theta[hyperparsel]
      prec.marginal   <- exp(theta[1L])
      rho.intern      <- theta[2L]
      rho             <- 2.0 * exp(rho.intern) / (1.0 + exp(rho.intern)) - 1.0
      prec.innovation <- prec.marginal/(1-rho^2)

      D <- bandSparse(N, N,
                      (-1):1,
                      list(rep(-rho, N-1),
                           rep(1, N),
                           rep(0, N-1)),
                      repr = "T")

      D[1,1] = sqrt(1-rho^2)

      return(sqrt(prec.innovation)*D)
    }
  }


  else if(model == "iid"){  #case when scale already present is not implemented
    D = as(Diagonal(N, 1),"dgCMatrix")
    Dfunc <- function(theta){
      prec <- exp(theta[hyperparsel])
      return(sqrt(prec)*D)
    }
  }


}


ar1.get.prior <- function(fit, name){

  components <- names(fit$summary.random)  #name of each component
  k <- which.max(name == components)       #find index for component name

  theta1 <- fit$all.hyper$random[[k]]$hyper$theta1
  theta2 <- fit$all.hyper$random[[k]]$hyper$theta2

  if( theta1$prior != "pc.prec") stop("ngvb currently only allows AR1 model with PC prior pc.prec")
  if( theta2$prior != "pc.cor0") stop("ngvb currently only allows AR1 model with PC prior pc.cor0")

  U1     <- theta1$param[1]
  alpha1 <- theta1$param[2]
  U2     <- theta2$param[1]
  alpha2 <- theta2$param[2]

  return(list(U1 = U1, alpha1 = alpha1, U2 = U2, alpha2 = alpha2))
}

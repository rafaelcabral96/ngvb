## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(ngvb)

## ---- fig.height=5, fig.width=15----------------------------------------------
plot(jumpts)

## -----------------------------------------------------------------------------
prec.prior <- list(prec = list(prior = "loggamma", param = c(1, 0.1)))

## -----------------------------------------------------------------------------
prior.list <- list(theta1 = list(prior = "pc.prec", param = c(1,0.1)),
                   theta2 = list(prior = "pc.cor0", param = c(0.5,0.5)))

## -----------------------------------------------------------------------------
formula <- y ~ -1 + f(x,  model = "ar1", hyper = prior.list) 
LGM     <- inla(formula, 
                data = jumpts,
                control.family = list(hyper = prec.prior))

## -----------------------------------------------------------------------------
LnGM <- ngvb(fit = LGM, selection = list(x=1:100))

## -----------------------------------------------------------------------------
#?`summary,ngvb.list-method`
#?`print,ngvb.list-method`
#?`plot,ngvb.list,missing-method`
#?`fitted,ngvb.list-method`
#?`simulate,ngvb.list-method`

## ---- fig.height=5, fig.width=15----------------------------------------------
plot(LnGM)

## -----------------------------------------------------------------------------
LGM$mlik[2]
LnGM@LGM$mlik[2]

## ---- fig.height=5, fig.width=15----------------------------------------------
plot(Orthodont.plot, layout = c(16,2))

## -----------------------------------------------------------------------------
summary(Orthodont)

## -----------------------------------------------------------------------------
formula <- value ~ 1 + Female + time + tF + f(subject, model = "iid") + f(subject2, time, model = "iid")

LGM <- inla(formula,
            data = Orthodont)

## ---- fig.height=5, fig.width=15----------------------------------------------
LnGM <- ngvb(LGM, selection = list(subject = 1:27, subject2 = 1:27))

plot(LnGM)

## -----------------------------------------------------------------------------
LGM$mlik[2]
LnGM@LGM$mlik[2]

## ----message=FALSE, warning=FALSE, fig.height=5, fig.width=15-----------------
#install.packages("spdep")
#install.packages("rgdal")
library(spdep)                # Columbus dataset
library(rgdal)                # Read polygon data

data(columbus)
data    <- columbus[,c("CRIME","HOVAL","INC")]       # data
N    <- nrow(data)                                   # number of counties
data$s  <- 1:N 
map  <- readOGR(system.file("shapes/columbus.shp", package="spData")[1]) #shape file containing the polygons
plot(map)


nb_q <- poly2nb(map)                                   # Construct neighbours list from polygon list
nb_W <- nb2listw(nb_q, style="B", zero.policy=TRUE)    # Spatial weights for neighbours lists
W    <- as.matrix(as(nb_W, "sparseMatrix"))            # Adjacency matrix W
W    <- diag(1/rowSums(W))%*%W                         # Row standardize adjacency matrix  
eigenv    <- eigen(W)$values                           # Eigenvalues of W. We need them to compute the log-determinant of the precision matrix

## -----------------------------------------------------------------------------
'inla.rgeneric.sar.model' <- function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL) {
  
  envir = parent.env(environment())
  
  #Internal scale is log of precision and log of kappa
  interpret.theta <- function() {
    return(list(prec = exp(theta[1L]),
                rho = exp(theta[2L])/(1+exp(theta[2L]))))
  }
  
  
  graph <- function(){
    return(Q())
  }
  
  Q <- function() { 
    param = interpret.theta()
    rho <- param$rho
    prec  <- param$prec
    
    D = Diagonal(n,rep(1,n)) - rho*W
    
    return(prec*t(D)%*%Diagonal(n,1/V)%*%D)
  }
  
  mu <- function()
  {
    return(numeric(0))
  }
  
  log.norm.const <- function() {
    
    param = interpret.theta()
    rho <- param$rho
    prec  <- param$prec
    
    log.det.D <- sum(log(1-rho*eigenv))
    
    res <- -0.5*n*log(2*pi) + 0.5*n*log(prec) + log.det.D -0.5*sum(log(V))
    return (res)  }
  
  

  log.prior <- function() {
    
    param = interpret.theta()
    rho <- param$rho
    prec  <- param$prec
    
    tau_prec = 1/sqrt(5)  
    
    prior.theta1 <- log( (2*sqrt(tau_prec/(2*pi)))*exp(-0.5*tau_prec*prec^2) ) + theta[1L]
    prior.theta2 <- log(exp(theta[2L])/(1+exp(theta[2L]))^2)
    
    res <-  prior.theta1 + prior.theta2
    
    return(res)
  }
  
  
  initial <- function() {
    return(c(0,0))
  }
  
  
  quit <- function() {
    return(invisible())
  }
  
  if (!length(theta)) theta = initial()
  res <- do.call(match.arg(cmd), args = list())
  return(res)
}

## -----------------------------------------------------------------------------
inla.fit.V <- function(V){
  
  model = inla.rgeneric.define(inla.rgeneric.sar.model, n = N, V = V[[1]], W = W, eigenv = eigenv)

  formula <- CRIME ~ 1 + HOVAL + INC  + f(s, model = model)

  fit <- inla(formula, 
              data = data,
              control.compute = list(config = TRUE))
  
  return(fit)
}

## -----------------------------------------------------------------------------
###Dfuncsel, h
selection <- list(s=1:N)
h         <- list(h=rep(1,N))

func1    <- function(theta){
  prec   <- exp(theta[2L]) 
  rho    <- exp(theta[3L])/(1+exp(theta[3L]))
  return(sqrt(prec)*(Diagonal(N,rep(1,N)) - rho*W))
}
Dfunc <- list(func1)

inla.fit.V <- function(V,theta=NULL){
  
  model = inla.rgeneric.define(inla.rgeneric.sar.model, n=N[[1]], V=V[[1]], W = W, eigenv = eigenv)

  formula <- CRIME ~ 1 + HOVAL + INC  + f(s, model = model)

  control.mode = NULL
  if(!is.null(theta)){
    control.mode = list(theta = theta, fixed = TRUE)
  }
  
  fit <- inla(formula, 
               control.compute = list(config = TRUE, cpo = TRUE),
               data = data, 
               control.mode = control.mode)
  
  return(fit)
}

manual.configs <- list(h = h, Dfunc = Dfunc, inla.fit.V = inla.fit.V)


## ---- fig.height=5, fig.width=15----------------------------------------------
LGM <- inla.fit.V(list(rep(1,N)))

LGM$summary.hyperpar
LGM$mlik[2]

#posterior marginals of the precision tau_x and  autocorrelation parameter rho
plot(inla.tmarginal(function(x)  exp(x), LGM$marginals.hyperpar$`Theta2 for s`), type = 'l', col = 'red', main = "precision")
plot(inla.tmarginal(function(x)  exp(x)/(1+exp(x)), LGM$marginals.hyperpar$`Theta2 for s`), type = 'l', col = 'red', main = "rho")

## -----------------------------------------------------------------------------
D1  <- function(theta){
  prec   <- exp(theta[2L]) 
  rho    <- exp(theta[3L])/(1+exp(theta[3L]))
  return(sqrt(prec)*(Diagonal(N,rep(1,N)) - rho*W))
}
Dfunc <- list(D1)

h         <- list(rep(1,N))

manual.configs <- list(inla.fit.V = inla.fit.V, Dfunc = Dfunc, h = h)

## ---- fig.height=5, fig.width=15----------------------------------------------
LnGM <- ngvb(manual.configs = manual.configs, selection = list(s=1:N), 
             iter = 10, d.sampling = TRUE, n.sampling = 1000)

## -----------------------------------------------------------------------------
plot(LnGM)
LnGM@LGM$mlik[2]

#posterior marginals of the precision tau_x and  autocorrelation parameter rho
plot(inla.tmarginal(function(x)  exp(x), LnGM@LGM$marginals.hyperpar$`Theta2 for s`), type = 'l', col = 'red', main = "precision")
plot(inla.tmarginal(function(x)  exp(x)/(1+exp(x)), LnGM@LGM$marginals.hyperpar$`Theta2 for s`), type = 'l', col = 'red', main = "rho")


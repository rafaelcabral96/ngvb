
##Returns prior in Log scale----------
eta.prior.f <- function(d,h,theta,N){
  f <- function(eta){
    logf <- numeric(length(eta))
    for(i in 1:length(eta)){
      value = sqrt(eta[i]^(-1)*(d+h^2/eta[i]))
      logf[i] <- sum(h)/eta[i] - theta*eta[i] -N/2*log(eta[i]) + sum(log(besselK(value,-1,expon.scaled = TRUE))-value) - 0.5*sum(log(d*eta[i] + h^2))
    }
    return(logf)
  }
  return(f)
}


sampler.inverseCDF <- function(logpdf, supp.min, supp.max, supp.points = 100000, n.samples = 1){

  #Find mode of pdf
  mode <- optimize(logpdf, interval = c(supp.min,supp.max), maximum = T)$maximum
  pdf <- function(x) exp(logpdf(x) - logpdf(mode)) #subtract by logpdf at mode to stabilize -> avoid overflow

  #Find region (pdfmin, pdfmax) that contains most of the mass
  imin = 1
  pdfmode <- pdf(mode)
  pdfmin  <- pdfmode
  while(pdfmin/pdfmode>1e-7){
    imin <- imin +1
    pointmin  <- mode/1.3^imin
    pdfmin <- pdf(pointmin)
  }

  imax = 1
  pdfmax  <- pdfmode
  while(pdfmax/pdfmode>1e-7){
    imax <- imax +1
    pointmax  <- mode*1.3^imax
    pdfmax <- pdf(pointmax)
  }


  #Do inverse cdf approximation on that region
  x    <- seq(pointmin, pointmax, length.out = supp.points)

  dens <- pdf(x)
  cdens <- cumsum(dens)
  cdens <- cdens/cdens[length(cdens)]

  duplicated <- cdens > 0.999999
  x          <- x[!duplicated]
  cdens      <- cdens[!duplicated]

  #TAKE N.SAMPLES FROM UNIFORM DISTRIBUTION
  r <- runif(n.samples)
  sample <- spline(x = cdens, y = x, xout = r)$y


  return( sample )
}


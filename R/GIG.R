rGIG <- function(n.sample, p, a, b){
  return(ngme2::rgig(p = rep(p, n.sample),
                    a = rep(a, n.sample),
                    b = rep(b, n.sample)))
}


mGIG <- function(p,a,b,order=1, n.samples = 1000){
  samples <- rGIG(n.samples,p,a,b)
  return(sum(samples^order)/n.samples)
}

#First moment of GIG distribution (p,a,b can be vectors)
GIGM1 <- function(p,a,b){
  return( sqrt(b/a)*besselK(sqrt(a*b),p+1)/besselK(sqrt(a*b),p) ) #
}

#First moment of GIG distribution (p,a,b can be vectors)
GIGM2 <- function(p,a,b){
  return( (b/a)*besselK(sqrt(a*b),p+2)/besselK(sqrt(a*b),p) ) #
}

#Moment of order -1 of GIG distribution (p,a,b can be vectors)
GIGMm1 <- function(p,a,b){
  return( sqrt(a/b)*besselK(sqrt(a*b),p+1)/besselK(sqrt(a*b),p) - 2*p/b)
}

#first mode of GIG distribution
GIGmode <- function(p,a,b){
  return( ((p-1)+sqrt((p-1)^2+a*b))/a )
}

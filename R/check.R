check.so.exists <- function(){
  if(!file.exists('cfiles/cgeneric-ngvb.so')) {

    folder = "cfiles/"
    ## generate the so
    system(paste0('gcc -c -fpic -I', folder, ' ',
                  folder, 'cgeneric-ngvb.c -o ', folder,'cgeneric-ngvb.o'))
    system(paste0('gcc -shared ', folder,'cgeneric-ngvb.o -o ', folder, 'cgeneric-ngvb.so'))
  }
}


change.args.controls <- function(r, fast){

  r$.args$control.compute$config <- TRUE
  r$.args$inla.mode <- "experimental"

  if(fast){
    r$.args$control.compute$return.marginals <- FALSE
    r$.args$control.compute$return.marginals.predictor <- FALSE
    r$.args$control.compute$dic <- FALSE
    r$.args$control.compute$cpo <- FALSE
    r$.args$control.compute$po <- FALSE
    r$.args$control.compute$waic <- FALSE
    r$.args$control.compute$graph <- FALSE
    r$.args$control.compute$hyperpar <- FALSE
    r$.args$control.compute$q <- FALSE

    r$.args$control.fixed$correlation.matrix <- FALSE

    r$.args$control.inla$int.strategy <- "eb"
    #r$.args$control.inla$use.directions <- r$misc$opt.directions

    #r$.args$control.mode$result <- NULL
    #r$.args$control.mode$restart <- FALSE
    #r$.args$control.mode$theta <- r$mode$theta
    #r$.args$control.mode$x <- r$mode$x
    #r$.args$control.mode$fixed <- TRUE

    r$.args$lincomb <- NULL
  }
  return(r)
}


ngvb.summary.mixing <- function(fit, configs){
  method    <- configs$method
  ncomp     <- configs$ncomp
  compnames <- names(configs$selection)
  h         <- configs$h

  if(method == "Gibbs"){
    V       <- fit$ngvb$V
    summary <-  data.frame(component = rep(compnames, lengths(V)),
                           index = unlist(lapply(lengths(V), function(x) 1:x)),
                           V = unlist(V))

  }
  else if(is.null(fit$ngvb$d)){ #first iteration
    V       <- fit$ngvb$V
    summary <- data.frame()
    for(k in 1:ncomp){
      for(i in 1:length(V[[k]])){
        summary <- rbind(summary, c(compnames[k], i, V[[k]][i]) )
      }
    }
    colnames(summary) <- c("component", "V index", "mean")
    cols     <- 2:3
    summary[cols] <- lapply(summary[cols], as.numeric)

  }
  else if(method == "SVI"){

    summary <- data.frame()
    j = 1
    for(k in 1:ncomp){
      list <- fit$ngvb$post.V[[k]]
      p    <- list$p_V
      a    <- list$a_V
      b    <- list$b_V
      Nk   <- length(b)

      mean <- GIGM1(p,a,b)
      sd   <- sqrt(GIGM2(p,a,b) - mean^2)
      mode <- GIGmode(p,a,b)


      data.k <- data.frame(rep(compnames[k],Nk), 1:Nk, mean, sd,
                           rep(NA,Nk), rep(NA,Nk), rep(NA,Nk),
                           mode, rep(p,Nk), rep(a,Nk), b)
      summary <- rbind(summary, data.k)

      #for(i in 1:length(b)){
      #  mean  = mGIG(p,a,b[i],order = 1)
      #  sd    = sqrt(mGIG(p,a,b[i],order = 2) - mean^2)
      #  q025  = qGIG(0.025,p,a,b[i])
      #  q50   = qGIG(0.50, p,a,b[i])
      #  q975  = qGIG(0.975,p,a,b[i])
      #  mode  = GIGmode(p,a,b[i])
      #  summary <- rbind(summary,c(compnames[k], i, mean, sd, q025, q50, q975, mode, p, a, b[i]))
      #  j = j + 1
      #}
    }

    colnames(summary) <- c("component", "V index", "mean", "sd", "0.025quant", "0.5quant",
                           "0.975quant", "mode", "GIG.p", "GIG.a", "GIG.b")
    cols     <- 2:11
    summary[cols] <- lapply(summary[cols], as.numeric)

  }
  else if(method == "SCVI"){
    summary <- data.frame()

    for(k in 1:ncomp){
      Vm    = fit$ngvb$SCVI.samples[[k]]
      mean  = apply(Vm, 2, mean)
      sd    = apply(Vm, 2, sd)
      q025  = apply(Vm, 2, quantile, 0.025)
      q50   = apply(Vm, 2, quantile, 0.5)
      q975  = apply(Vm, 2, quantile, 0.975)
      summary = rbind(summary, cbind(compnames[k], 1:ncol(Vm), mean, sd, q025, q50, q975))
    }
    colnames(summary) <- c("component", "V index", "mean", "sd", "0.025quant", "0.5quant", "0.975quant")
    cols     <- 2:7
    summary[cols] <- lapply(summary[cols], as.numeric)
  }
  else{ stop("Error: method should be either SVI or SCVI")}

  return(summary)
}

ngvb.summary.ng <- function(fit, configs){
  method   <- configs$method
  ncomp    <- configs$ncomp

  if(method == "Gibbs"){
    summary <- data.frame(eta = fit$ngvb$eta)
  }
  else if(is.null(fit$ngvb$d)){
    summary <- matrix(NA, nrow = ncomp, ncol = 1)
    colnames(summary) <- c("mean")
    for(k in 1:ncomp){
      mean        <- fit$ngvb$eta[k]
      summary[k,] <- c(mean)
    }
  }
  else if(method == "SVI"){
    summary <- matrix(NA, nrow = ncomp, ncol = 9)
    colnames(summary) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode", "GIG.p", "GIG.a", "GIG.b")
    for(k in 1:ncomp){
      list  <- fit$ngvb$post.eta[[k]]
      p     <- list$p_eta
      a     <- list$a_eta
      b     <- list$b_eta
      samples <- rGIG(10000, p, a, b)
      mean  <- mean(samples)
      sd    <- sd(samples)
      q025  <- quantile(samples, 0.025)
      q50   <- quantile(samples, 0.50)
      q975  <- quantile(samples, 0.975)
      mode  <- GIGmode(p,a,b)
      summary[k,] <- c(mean, sd, q025, q50, q975, mode, p, a, b)
    }

  }
  else if(method == "SCVI"){
    summary <- matrix(NA, nrow = ncomp, ncol = 5)
    colnames(summary) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
    for(k in 1:ncomp){
      samples     <- fit$ngvb$eta.samples[[k]]
      mean        <- mean(samples)
      sd          <- sd(samples)
      q025        <- quantile(samples, 0.025)
      q50         <- quantile(samples, 0.50)
      q975        <- quantile(samples, 0.975)
      summary[k,] <- c(mean, sd, q025, q50, q975)
    }

  }
  else{ stop("Error: method should be either SCVI or SVI") }

  summary           <- as.data.frame(summary)
  rownames(summary) <- sapply(names(configs$selection), function(x) paste0("Non Gaussianity parameter for ", x) )

  return(summary)
}


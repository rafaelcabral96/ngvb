`ngvb` is an R package that performs inference of latent non-Gaussian models using variational Bayes and Laplace approximations. Several use cases can be found in the online [vignette](https://htmlpreview.github.io/?https://github.com/rafaelcabral96/ngvb/blob/master/vignettes/ngvb.html).


It requires the packages `inla` and `ngme` which can be installed by:

```
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
remotes::install_github("davidbolin/ngme", ref = "devel")
```

The `ngvb` package can be installed using the command:

```
devtools::install_github("rafaelcabral96/ngvb", build_vignettes = TRUE)
```

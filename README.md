`ngvb` is an R package that performs inference of latent non-Gaussian models using variational Bayes and Laplace approximations. Several use cases can be found in the online [vignette](http://htmlpreview.github.io/?https://raw.githubusercontent.com/rafaelcabral96/ngvb/master/vignettes/ngvb.html?token=GHSAT0AAAAAABYQK4RHOAAJ7VU3AJR2OTXWY3I3ZOQ).


It requires the packages `INLA` and `ngme` which can be installed by:

```
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
remotes::install_github("davidbolin/ngme", ref = "devel")
```

The `ngvb` package can be installed using the command:

```
devtools::install_github("rafaelcabral96/ngvb", build_vignettes = TRUE)
```

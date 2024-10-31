
# About

`ngvb` is an R package that performs inference of latent non-Gaussian models using variational Bayes and Laplace approximations. Read the  [Get started to the ngvb package](https://rafaelcabral96.github.io/ngvb/articles/ngvb.html) for an introduction to the package with several examples.

# Installation instructions

It requires the packages `INLA` and `ngme` which can be installed by:

```
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
remotes::install_github("davidbolin/ngme2", ref = "devel")
```

For more details on installation, visit the [ngme2 webpage](https://davidbolin.github.io/ngme2/) and the [R-INLA webpage](https://www.r-inla.org/).

The `ngvb` package can be installed using the command:

```
devtools::install_github("rafaelcabral96/ngvb", build_vignettes = TRUE)
```

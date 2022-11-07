`ngvb` is an R package that performs inference of latent non-Gaussian models using variational Bayes and Laplace approximations. 

It requires the packages `inla` and `ngme` which can be installed by:

```
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
remotes::install_github("davidbolin/ngme", ref = "devel")
```

The development version of the `ngvb` package can be installed using the command

```
devtools::install_github("rafaelcabral96/ngvb")
```

After installing, run `devtools::build_vignettes("ngvb")` and `vignette("ngvb")` to see several use cases.

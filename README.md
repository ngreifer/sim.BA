
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simBA

<!-- badges: start -->
<!-- badges: end -->

Unmeasured confounding is often raised as a source of potential bias
when evaluating non-randomized study protocols, but evaluating such
concerns during their design remains challenging. We propose a flexible
methodology based on individual level simulations that can allow
researchers to characterize the bias arising from unmeasured confounding
with a specified but modifiable structure during the study design.

`simBA` allows user to conduct a simulation-based quantitative bias
analysis using covariate structures generated with individual-level data
to characterize the bias arising from unmeasured confounding. Users can
specify their desired data generating mechanisms to simulate data and
quantitatively summarize findings in an end-to-end application using
this package. See `vignette("simBA")` for details.

## Installation

You can install the development version of `simBA` from GitLab with:

``` r
# install.packages("remotes")
remotes::install_gitlab("rjd48/simBA", host = "gitlab-scm.partners.org")
```

You can install the published version from CRAN with:

``` r
install.packages("simBA")
```

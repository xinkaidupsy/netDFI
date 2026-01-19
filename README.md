
# Dynamic Fit Index Cutoffs for Network Models
<!-- badges: start -->
  [![R-CMD-check](https://github.com/xinkaidupsy/netDFI/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/xinkaidupsy/netDFI/actions/workflows/R-CMD-check.yaml)
  [![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
  [![CRAN status](https://www.r-pkg.org/badges/version/netDFI)](https://CRAN.R-project.org/package=netDFI)
  [![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/netDFI)](https://CRAN.R-project.org/package=netDFI)
<!-- badges: end -->

The goal of `netDFI` is to compute dynamic fit index (DFI) cutoffs for network models. 

The reference of DFI for Gaussian graphical models: https://osf.io/preprints/psyarxiv/5wj2y_v1

## Installation

You can install the development version of netDFI from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xinkaidupsy/netDFI")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(psych)
library(psychonetrics)
library(dplyr)
library(netDFI)

# get the big five inventory data from psych
data(bfi)

# estimate ggm 
bfi_mod <- ggm(bfi) %>% prune %>% runmodel     

# obtain the partial correlation matrix
bfi_net <- getmatrix(bfi_mod, "omega")

# allow future_apply to use more memory
options(future.globals.maxSize = 2 * 1024^3)

# run dfi
dfi_bfi <- dfi_ggm(
  bfi_net, 
  ncores = parallel::detectCores(), 
  power = 0.80, 
  iter = 200,
  n_misspec = 2
)
dfi_bfi

# plot results
p <- plot(dfi_bfi)
p[[1]]

```


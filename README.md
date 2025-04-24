
# netDFI

<!-- badges: start -->
<!-- badges: end -->

The goal of netDFI is to ...

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

data(bfi)

bfi_mod <- ggm(bfi) %>% prune %>% runmodel     

bfi_net <- getmatrix(bfi_mod, "omega")

dfi_bfi <- dfi_ggm(bfi_net)
```


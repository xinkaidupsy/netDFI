
# netDFI

<!-- badges: start -->
<!-- badges: end -->

The goal of netDFI is to compute dynamic fit index cutoffs for network models. 

[Gaussian graphical models publication]{https://osf.io/preprints/psyarxiv/5wj2y_v1}

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
options(future.globals.maxSize = 1024 * 1024^2)

# run dfi
dfi_bfi <- dfi_ggm(bfi_net, ncores = 8, power = 0.95, iter = 500)

```


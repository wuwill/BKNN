
# BKNN

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of BKNN is to implement the Baysian k-Nearest Neighbhors algorithm in R.

## Installation

You can install the development version of BKNN like so:

``` r
devtools::install_github("wuwill/BKNN")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(BKNN)
dat <- matrix() # a data set with 10 samples
bnn_res <- calc_bnn(dat, target_sample = colnames(dat)[1:3])
```

## Reference
Nuti G. An efficient algorithm for bayesian nearest neighbours. Methodology and Computing in Applied Probability. 2019 Dec;21(4):1251-8.

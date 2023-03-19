
# BKNN

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of BKNN is to implement the Bayesian k-Nearest Neighbors algorithm in R.

## Installation

You can install the development version of BKNN like so:

``` r
devtools::install_github("wuwill/BKNN")
```

## Example

``` r
library(BKNN)
dat <- cbind(matrix(rnorm(30, 0), ncol = 5),
        matrix(rnorm(30, 2), ncol = 5))
sample_info <- data.frame(Name = paste0("I", 1:10), gender = rep(c("MALE", "FEMALE"), each = 5))
bnn_res <- calc_bnn(dat, sample_info, query = sample_info$Name, target = sample_info$Name, sample_col = "Name", status_col = "gender", class1 = "MALE")
```

## Reference
Nuti G. An efficient algorithm for Bayesian nearest neighbours. Methodology and Computing in Applied Probability. 2019 Dec;21(4):1251-8.

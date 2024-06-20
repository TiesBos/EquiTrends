
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EquiTrends

<!-- badges: start -->

[![R-CMD-check](https://github.com/TiesBos/EquiTrends/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TiesBos/EquiTrends/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Testing for parallel trends is crucial in the Difference-in-Difference
framework. The goal of EquiTrends is to provide a set of functions to
test for equivalence of pre-trends in difference-in-differences
estimation. The procedures follow the work of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121)). The package
provides functions to test for equivalence of pre-trends in
difference-in-differences estimation using the placebo coefficient
estimates, used to compare the trends in the pre-treatment period to
some base period (generally the final period in the pre-treatment
period). The package contains the functions to perform the testing
procedure surrounding the maximum placebo coefficient (see equation
(3.1) of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121))), to perform the
testing procedure surrounding the mean placebo coefficient (see equation
(3.2) of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121))) and to perform
the testing procedure surrounding the root mean squared placebo
coefficient (see equation (3.3) and (3.4) of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121))). Furthermore,
the package contains the function to simulate a paneldataset for such
testing purposes.

## Installation

You can install the development version of EquiTrends from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("TiesBos/EquiTrends")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
#library(EquiTrends)
## basic example code
```


<!-- README.md is generated from README.Rmd. Please edit that file -->

# drugdevelopR: Utility-based optimal phase II/III drug development planning.

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/Sterniii3/drugdevelopR.svg?branch=master)](https://travis-ci.com/Sterniii3/drugdevelopR)
[![CRAN
status](https://www.r-pkg.org/badges/version/drugdevelopR)](https://CRAN.R-project.org/package=drugdevelopR)
[![R-CMD-check](https://github.com/Sterniii3/drugdevelopR/workflows/R-CMD-check/badge.svg)](https://github.com/Sterniii3/drugdevelopR/actions)
<!-- badges: end -->

The drugdevelopR package enables you to plan phase II/III drug
development programs with optimal sample size allocation and go/no-go
decision rules. The assumed true treatment effects can be fixed or
modelled by a prior distribution. The corresponding [R Shiny
application](https://web.imbi.uni-heidelberg.de/drugdevelopR/) has a
graphic user interface for the package and thus makes it accessible for
users without prior knowledge of R. Fast computing is made possible by
parallel programming. The theoretical foundations for this package were
laid in the dissertation “Integrated Planning of Pilot and Subsequent
Confirmatory Study in Clinical Research – Finding Optimal Designs in a
Utility-Based Framework” by Stella Erdmann at the Institute of Medical
Biometry at the University of Heidelberg.

On our webpage, we supply [full
documentation](https://sterniii3.github.io/drugdevelopR/references/index.html)
of all functions as well as a [tutorial for getting
started](https://sterniii3.github.io/drugdevelopR/vignettes/this-link-is-not-yet-working)
with drugdevelopR.

## Installation

Install the development version of the package directly from
[GitHub](https://github.com/Sterniii3/drugdevelopR/) using the following
code:

``` r
if(!require(devtools)) { install.packages("devtools"); require(devtools)} 
devtools::install_github("Sterniii3/drugdevelopR")
```

and access the drugdevelopR App via
<https://web.imbi.uni-heidelberg.de/drugdevelopR/>.

## Usage

Here is a basic example for applying drugdevelopR to a drug development
program with a normally distributed outcome:

``` r
library(drugdevelopR)
#> Loading required package: mvtnorm
#> Loading required package: doParallel
#> Loading required package: foreach
#> Loading required package: iterators
#> Loading required package: parallel
#> Loading required package: msm
#> Loading required package: cubature
# TODO: Fill example here.
```

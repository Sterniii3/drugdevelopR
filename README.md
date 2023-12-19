
<!-- README.md is generated from README.Rmd. Please edit that file -->

# drugdevelopR: Utility-based optimal phase II/III drug development planning.

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/drugdevelopR)](https://CRAN.R-project.org/package=drugdevelopR)
[![R-CMD-check](https://github.com/Sterniii3/drugdevelopR/workflows/R-CMD-check/badge.svg)](https://github.com/Sterniii3/drugdevelopR/actions)
[![Codecov test
coverage](https://codecov.io/gh/Sterniii3/drugdevelopR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Sterniii3/drugdevelopR?branch=master)
<!-- badges: end -->

The drugdevelopR package enables you to plan phase II/III drug
development programs with optimal sample size allocation and go/no-go
decision rules. The assumed true treatment effects can be fixed or
modelled by a prior distribution. The corresponding [R Shiny
application](https://web.imbi.uni-heidelberg.de/drugdevelopR/) has a
graphic user interface and thus makes it accessible for
users without prior knowledge of R. Fast computing is made possible by
parallel programming. Theoretical foundations for this package were laid
in the dissertation “Integrated Planning of Pilot and Subsequent
Confirmatory Study in Clinical Research – Finding Optimal Designs in a
Utility-Based Framework” by Stella Erdmann at the Institute of Medical
Biometry at the University of Heidelberg.

On the package webpage, we supply [full
documentation](https://sterniii3.github.io/drugdevelopR/reference/index.html)
of all functions as well as a [tutorial for getting
started](https://sterniii3.github.io/drugdevelopR/articles/Introduction-to-drugdevelopR.html)
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
#> Lade nötiges Paket: doParallel
#> Lade nötiges Paket: foreach
#> Lade nötiges Paket: iterators
#> Lade nötiges Paket: parallel
#> Lade nötiges Paket: progressr
# Optimize
optimal_normal(Delta1 = 0.625, Delta2 = 0.8, fixed = FALSE, # treatment effect
               n2min = 20, n2max = 400, # sample size region
               stepn2 = 4, # sample size step size
               kappamin = 0.02, kappamax = 0.2, # threshold region
               stepkappa = 0.02, # threshold step size
               c2 = 0.675, c3 = 0.72, # maximal total trial costs
               c02 = 15, c03 = 20, # maximal per-patient costs
               b1 = 3000, b2 = 8000, b3 = 10000, # gains for patients
               alpha = 0.025, # one-sided significance level
               beta = 0.1, # 1 - power
               w = 0.6, in1 = 300, in2 = 600, # weight and amount of information
               a = 0.25, b = 0.75) # truncation values
#> Optimization result:
#>  Utility: 3073.85
#>  Sample size:
#>    phase II: 84, phase III: 194, total: 278
#>  Probability to go to phase III: 0.99
#>  Total cost:
#>    phase II: 72, phase III: 160, cost constraint: Inf
#>  Fixed cost:
#>    phase II: 15, phase III: 20
#>  Variable cost per patient:
#>    phase II: 0.675, phase III: 0.72
#>  Effect size categories (expected gains):
#>   small: 0 (3000), medium: 0.5 (8000), large: 0.8 (10000)
#>  Success probability: 0.85
#>  Success probability by effect size:
#>    small: 0.69, medium: 0.15, large: 0.01
#>  Significance level: 0.025
#>  Targeted power: 0.9
#>  Decision rule threshold: 0.06 [Kappa] 
#>  Parameters of the prior distribution: 
#>    Delta1: 0.625, Delta2: 0.8, in1: 300, in2: 600,
#>    a: 0.25, b: 0.75, w: 0.6
#>  Treatment effect offset between phase II and III: 0 [gamma]
```

## drugdevelopR functions

The drugdevelopR package provides the functions

- [`optimal_tte`](https://sterniii3.github.io/drugdevelopR/reference/optimal_tte.html),
- [`optimal_binary`](https://sterniii3.github.io/drugdevelopR/reference/optimal_binary.html),
  and
- [`optimal_normal`](https://sterniii3.github.io/drugdevelopR/reference/optimal_normal.html)

to plan optimal phase II/III drug development programs with

- time-to-event (treatment effect measured by hazard ratio, *HR*),
- binary (treatment effect measured by risk ratio, *RR*), or
- normally distributed (treatment effect measured by standardized
  difference in means, *Delta*)

endpoints, where the treatment effect is modelled on a [prior
distribution](https://web.imbi.uni-heidelberg.de/prior/). Optimal phase
II/III drug development planning with fixed treatment effects can be
done with the help of the R Shiny application
[basic](https://web.imbi.uni-heidelberg.de/basic/).

Extensions to the basic setting are:

- optimal planning of programs including methods for discounting of
  phase II results (function:
  [optimal_bias](https://sterniii3.github.io/drugdevelopR/reference/optimal_bias.html),
  App: [bias](https://web.imbi.uni-heidelberg.de/bias/)),
- optimal planning of programs with several phase III trials (function:
  [optimal_multitrial](https://sterniii3.github.io/drugdevelopR/reference/optimal_multitrial.html),
  App: [multitrial](https://web.imbi.uni-heidelberg.de/multitrial/)) and
- optimal planning of programs with multiple arms (function:
  [optimal_multiarm](https://sterniii3.github.io/drugdevelopR/reference/optimal_multiarm.html),
  App: [multiarm](https://web.imbi.uni-heidelberg.de/multiarm/)).

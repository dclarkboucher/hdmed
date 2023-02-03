
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hdmed

<!-- badges: start -->
<!-- badges: end -->

Our package offers a suite of functions for performing mediation
analysis with high-dimensional mediators. Unlike methods for
single-mediator mediation analysis—which have been distributed by
packages such as “[psych](https://CRAN.R-project.org/package=psych),”
“[mediation](https://CRAN.R-project.org/package=mediation),”
“[medScan](https://CRAN.R-project.org/package=medScan)”—our package
focuses on settings whether there are many potential mediators that need
evaluating simultaneously, a topic which has recently become the focus
of prolific and exciting methodological work.

## Installation

You can install the development version of hdmed from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dclarkboucher/hdmed")
```

## Example

To demonstrate high-dimensional mediation analysis with an example, let
$\mathbf{M}$ be a set of $p$ variables, each a potential mediator in the
causal pathway between exposure $A$ and outcome $Y$, and let
$\mathbb{C}$ be a set of $q$ covariates, $q$ reasonably small. Then,
given data on individuals $i\in\{1,\dots,n\}$, with $n$ potentially but
not necessarily less than $p$, we can evaluate the mediating role of $M$
with the models

``` r
library(hdmed)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

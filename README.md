
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
$A$ be an exposure, $Y$ an outcome, $\mathbf{C}$ a set of $q$
covariates, and $\mathbf{M}$ a set of $p$ potential mediators in the
causal pathway between $A$ and $Y$. Then, supposing we have data on $n$
individuals, we can evaluate the mediating role of $\mathbf{M}$ with the
equations

$$
\begin{equation}
E[\mathbf{Y_i}|A_i,\mathbf{M}_i,\mathbf{C_i}] = \beta_aA_i+\mathbf{\beta_m}^T\mathbf{M_i} + \mathbf{\beta_c}^T\mathbf{C_i} 
\end{equation}
$$

and

$$
\begin{gather}
E[\mathbf{Y_i}|A_i,\mathbf{M}_i,\mathbf{C_i}] = \beta_aA_i+\mathbf{\beta_m}^T\mathbf{M_i} + \mathbf{\beta_c}^T\mathbf{C_i} 
\end{gather}
$$

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

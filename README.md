
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

You can install hdmed from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dclarkboucher/hdmed")
```

## Overview

To see how high-dimensional mediation analysis works mathematically, let
$A$ be an exposure, $Y$ be an outcome, $\mathbf{C}$ be a set of $q$
covariates, and $\mathbf{M}$ be a set of $p$ potential mediators in the
causal pathway between $A$ and $Y$. Then, supposing we have data on $n$
individuals, we can evaluate the mediating role of $\mathbf{M}$ with the
equations

$$
\begin{equation}
E[Y_i|A_i,\mathbf{M}_i,\mathbf{C_i}] = \beta_aA_i+\mathbf{\beta_m}^T\mathbf{M_i} + \mathbf{\beta_c}^T\mathbf{C_i} 
\end{equation}
$$

and

$$
\begin{equation}
E[\mathbf{M_i}|A_i,\mathbf{C_i}] =\mathbf{\alpha_a}A_i + \mathbf{\alpha_c}\mathbf{C_i}\text{,}
\end{equation}
$$

where the first equation is the **outcome model** and the second
equation is the **mediator model**. In the outcome model, our primary
estimands are $\beta_a$— the direct effect of the exposure on the
outcome independent of $\mathbf{M}$— and $\mathbf{\beta_m}$, a
$p$-vector of the association between each mediator and $Y$ given $A$
and $\mathbf{C}$. (Unlike many methods common to single-mediator
analysis, all the methods included in our package assume there is no
interaction effect between $\mathbf{M}$ and $A$ on $Y$.) Likewise, our
primary concern in the mediator model is $\mathbf{\alpha_a}$, which is a
$p$-vector of the conditional associations between each mediator and the
exposure given $\mathbf{C}$. As for the other coefficients,
$\mathbf{\beta_c}$ is a $q$-vector of the covariate-outcome effects, and
$\mathbf{\alpha_c}$ is a $p\times q$ matrix of covariate-mediator
associations.

Once the outcome and mediator models have been fitted, mediation
analysis can be performed by assessing their estimated coefficients. The
chief quantities of interest are:

1.  $\mathbf{\alpha_a}^T \mathbf{\beta_m}$, the **global mediation
    effect** of $A$ on $Y$ through $M$;

2.  $\beta_a$, the **direct effect** of $A$ on $Y$;

3.  $\mathbf{\alpha_a}^T \mathbf{\beta_m} + \beta_a$, the total affect
    of $A$ on $Y$; and

4.  $\frac{\mathbf{\alpha_a}^T \mathbf{\beta_m}}{\mathbf{\alpha_a}^T \mathbf{\beta_m}+\beta_a}$,
    the proportion of the total effect due to mediation (referred to as
    the *proportion mediated*.)

All the methods provided by our package can fit this model except for
HDMM (`mediate_hdmm`) and LVMA (`mediate_lvma`), which instead of the
typical model assumptions, assume the mediation between $A$ and $Y$ is
transmitted by unmeasured latent variables. (See the documentation of
those functions for more detail.) The other methods produce, at the very
least, estimates of the direct effect, global mediation effect, and
total effect, making them suitable for performing mediation analysis
with the standard assumptions. Moreover, in the case of BSLMM
(`mediate_bslmm`), HIMA (`mediate_hima`), HDMA (`mediate_hdma`), MedFix
(`mediate_medfix`), and Pathaway LASSO (`mediate_pathway_lasso`), we
also report estimates of the **mediation contributions**, which are the
contributions $(\mathbf{\alpha_a})_j(\mathbf{\beta_m})_j$ of each
mediator to $\mathbf{\alpha_a}^T \mathbf{\beta_m}$, $j$ from $1$ to $p$.
Though useful for identifying potentially important mediators, we stress
that these contributions *generally cannot be interpreted as causal
effects unless the mediators are* *independent conditional on* $A$ *and*
$\mathbf{C}$. Conditions for when $\mathbf{\alpha_a}^T \mathbf{\beta_m}$
and $\beta_a$ can be interpreted causally are laid out by Song et
al. (2019) (see `mediate_bslmm` for complete reference). Note also that,
as programmed, the methods HIMA (`mediate_hima`), HDMA (`mediate_hdma`),
MedFix (`mediate_medfix`), and BSLMM allow one to incorporate a small
number of covariates directly, as specified in the above pair of models,
whereas the other as programmed methods do not. If you are interested in
adjusting for covariates with a method that does not allow them to be
inputted to our mediation function directly, consider regressing those
covariates out of the outcome, mediators, and exposures in advance, when
doing so is appropriate. In addition, most functions in our package
assume that the outcome variable is continuous; however, HIMA and HDMA
have options for fitting a binary outcome model with a standard logistic
link.

## Example

<!-- badges: start -->
<!-- badges: end -->

The `med_dat` object provided by our package contains a simple toy
dataset for practicing high-dimensional mediation (though in this case,
we are using “high-dimensional” generously, as the dataset contains only
20 mediators to its 100 observations).

Let us take a look at the data. In `Y` we have the outcome, in `A` we
have the exposure, and in `M` we have a named matrix of mediators.

``` r
library(hdmed)

# Process data
Y <- med_dat$Y
M <- med_dat$M
A <- med_dat$A

str(M)
#>  num [1:100, 1:20] -0.491 1.339 -0.194 -0.218 -0.108 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:20] "m1" "m2" "m3" "m4" ...
```

Now we will perform mediation analysis. For a simple, fast mediation
method we will use the method “high-dimensional mediation analysis” by
Zhang et al. (2016), which we call “HIMA”. HIMA is a straightforward
method that fits the mediator models using ordinary least squares and
the outcome model using penalized regression with the minimax concave
penalty. We don’t have covariates to include, so to use the default
options, we input only `A`, `M`, and `Y`.

``` r
hima_out <- mediate_hima(A, M, Y)
```

Next let’s look at the mediation contributions, which are located in the
`contributions` table. In this case, the function only returned one
mediator, which happens if the others have an estimated contribution of
zero and do not contribute to the estimated global mediation effect.
Examining the table further, we see `alpha` as a shorthand for
$(\mathbf{\alpha_a})_j$, `beta` as a shorthand for
$(\mathbf{\beta_m})_j$, and `alpha_beta` as a shorthand for
$(\mathbf{\alpha_a})_j(\mathbf{\beta_m})_j$. Notice that `ab_pv` is the
$(\mathbf{\alpha_a})_j(\mathbf{\beta_m})_j$ p-value.

``` r
hima_out$contributions
#>   mediator      alpha   alpha_pv      beta      beta_pv alpha_beta      ab_pv
#> 1       m3 -0.2737383 0.01438887 0.6040214 1.364285e-07 -0.1653438 0.01438887
```

Finally, the estimated mediation effects are reported in `effects`
table, which includes the indirect effect (the global mediation effect),
the direct effect, and the total effect. In theory, one can use these
estimates to report the proportion mediated, as described above, but
since the proportion mediated is generally only useful when
$\mathbf{\alpha_a}^T \mathbf{\beta_m}$ and $beta_m$ have the same sign,
we will not do so here.

``` r
hima_out$effects
#>     effect    estimate
#> 1 indirect -0.16534380
#> 2   direct  0.01018211
#> 3    total -0.15516169
```

## Citation

This package serves as companion code for our paper, “Methods for
Mediation Analysis with High-Dimensional DNA Methylation Data: Possible
Choices and Comparison.” To give our work proper credit, use the
citation provided below:

> Dylan Clark-Boucher, Xiang Zhou, Jiacong Du, Yongmei Liu, Belinda L
> Needham, Jennifer A Smith, Bhramar \>Mukherjee (2023). Methods for
> Mediation Analysis with High-Dimensional DNA Methylation Data:
> Possible \>Choices and Comparison. medRxiv 2023.02.10.23285764

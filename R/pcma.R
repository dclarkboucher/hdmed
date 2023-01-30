#' Principal Component Mediation Analysis for High-dimensional Mediators
#'
#' @description \code{mediate_pcma} applies principal component mediation analysis
#' (Huang and Pan, 2013) to mediation settings in which the mediators are high-dimensional.
#'
#' @param A length \code{n} numeric vector containing exposure variable
#' @param M \code{n x p} numeric matrix of high-dimensional mediators.
#' @param Y length \code{n} numeric vector containing continuous outcome variable.
#' @param var_per a numeric variable with the desired proportion of variance
#' explained. Default is 0.8.
#' @param n_pc optional numeric variable with the desired number of PCs, in which case
#' \code{var_per} is ignored. Default is \code{NULL} and the number of PCs is
#' determined based on the desired proportion of variance explained.
#' @param sims number of Monte Carlo draws for nonparametric bootstrap or
#' quasi-Bayesian approximation (see [mediation::mediate()]).
#' Default is 1000.
#' @param boot_ci_type 	a character string indicating the type of bootstrap
#' confidence intervals for when \code{boot = TRUE}. If \code{"bca"},
#' bias-corrected and accelerated (BCa) confidence intervals will be estimated.
#' If \code{"perc"}, percentile confidence intervals will be estimated
#' (see [mediation::mediate()]). Default is "bca".
#' @param ci_level the desired confidence level. Default is 0.95.
#'
#' @details
#' Principal component mediation analysis (PCMA) is a method for estimating
#' mediation effects when the mediators are high-dimensional. The first step
#' is to compute the residuals of mediator models (\eqn{M|A}), then perform
#' PCA on those residuals to reduce them to a smaller number of mediators
#' that efficiently explain the residual variance. Then, since those mediators
#' are linearly independent conditional on A, one can trivially perform
#' single-mediator mediation analysis for each PC on its own, in this case
#' by using the [mediation::mediate()] function. The global mediation effect is estimated
#' by summing the mediation effects of the individual PCs.
#'
#'
#' @return A list containing:
#' \itemize{
#'     \item{loadings: }{a matrix of the PC loadings.}
#'     \item{pcs: }{a matrix of the PCs.}
#'     \item{var_explained: the cumulative proportion of variance explained by the PCs.}
#'     \item{contributions: }{a data frame containing the estimates, confidence
#'     intervals, and p-values of the mediation contributions.}
#'     \item{effects: }{a data frame containing the estimated direct, global
#'     mediation, and total effects}
#' }
#'
#'
#'
#' @import MASS
#' @import Matrix
#' @import mvtnorm
#' @import sandwich
#' @import mediation
#'
#' @references Huang, Y.-T. & Pan, W.-C. Hypothesis test of mediation effect in
#' causal mediation model with  high-dimensional continuous mediators.
#' Biometrics 72, 402–413 (2016).
#'
#'
#' @examples
#' data("med_dat")
#' A <- med_dat$A
#' M <- med_dat$M
#' Y <- med_dat$Y
#' # Fit PCMA with 3 principal components and print the effects. In practice one
#' # should choose n_pc (or var_per) and the number sims to be larger
#' out <- mediate_pcma(A, M, Y, n_pc = 3, sims = 10)
#' out$effects
#'
#'
#'
#' @export
#'

mediate_pcma <-
  function(A, M, Y, var_per = 0.8, n_pc = NULL, sims = 1000,
           boot_ci_type = "bca", ci_level = 0.95){

    # Create column names if absent
    p <- ncol(M)
    if (is.null(colnames(M))) colnames(M) <- paste0("m", 1:p)

    #Check A, M, Y
    if(is.data.frame(M)) M <- as.matrix(M)
    if(!is.numeric(A) | !is.vector(A)) stop("A must be numeric vector.")
    if(!is.numeric(M) | !is.matrix(M)) stop("M must be numeric matrix.")
    if(!is.numeric(Y) | !is.vector(Y)) stop("Y must be numeric vector.")
    if(is.null(colnames(M))){
      colnames(M) <- paste0("m",1:p)
    }

    # Check n_pc
    adaptive <- T
    if (!is.null(n_pc)){
      var_per <- NULL
      adaptive <- F
      message("Number of PCs provided; variance percentage ignored.")
    }

    # Obtain principal component loadings
    pca_out <-
      pca_loadings(A, M, adaptive = adaptive, var.per = var_per, n.pc = n_pc)

    # Organize loadings
    n_pc <- ncol(pca_out$U)
    loadings <- pca_out$U
    rownames(loadings) <- colnames(M)
    colnames(loadings) <- paste0("u",1:n_pc)

    # Compute PCs
    pcs <- M %*% loadings
    colnames(pcs) <- paste0("pc",1:n_pc)

    # Run marginal mediation on sparse PCs
    pcma_out <-
      mediate_multiple(
        A,
        pcs,
        Y,
        sims = sims,
        boot.ci.type = boot_ci_type,
        conf.level = ci_level
      )

    # Organize mediation contributions
    contributions <-
      with(
        pcma_out,
        data.frame(
          mediator = colnames(pcs),
          alpha = alpha[,1],
          beta = beta[,1],
          alpha_beta = IE[,1],
          ab_se = IE[,2],
          cl1 = IE[,4],
          cl2 = IE[,5],
          ab_pv = IE[,3]
        )
      )
    percentiles <- round(c((1 - ci_level) /2, 1 - (1 - ci_level) / 2) * 100,1)
    ci_names <-  paste0("ab_cl_", percentiles,"%")
    colnames(contributions)[6:7] <- ci_names

    # Organize effects
    effects <- matrix(NA, 3, 6)
    effects[, 1] <- c("indirect","direct","total")
    effects[2, 2:6] <- pcma_out$DE[1, c(1, 2, 4:5, 3)]
    effects[1, 2:6] <- pcma_out$IE.total[1, c(1, 2, 4:5, 3)]
    effects[3, 2:6] <- pcma_out$TE[1, c(1, 2, 4:5, 3)]
    ci_names <- paste0("cl_",percentiles,"%")
    colnames(effects) <- c("effect","estimate","se",ci_names,"pv")

    # Return output
    out <-
      list(
        loadings = as.data.frame(loadings),
        pcs = as.data.frame(pcs),
        var_explained = pca_out$var.pc,
        contributions = contributions,
        effects = as.data.frame(effects)
      )

    return(out)
  }

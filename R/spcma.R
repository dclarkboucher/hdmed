#' Sparse Principal Component Mediation Analysis for High-Dimensional Mediators
#'
#' @description \code{spcma} applies sparse principal component mediation
#' analysis to mediation settings in which the mediators are high-dimensional.
#'
#' @param A length \code{n} numeric vector containing exposure variable
#' @param M \code{n x p} numeric matrix of high-dimensional mediators.
#' @param Y length \code{n} numeric vector containing continuous outcome variable.
#' @param var_per a numeric variable with the desired proportion of variance
#' explained if \code{adaptive = TRUE}. Default is 0.8.
#' @param adaptive a logical variable for whether the number of PCs should be
#' determined based on the proportion of variance explained. Default is \code{TRUE}.
#' @param n_pc a numeric variable with the desired number of PCs for when
#' \code{adaptive = FALSE}. Default is \code{NULL}.
#' @param sims number of Monte Carlo draws for nonparametric bootstrap or
#' quasi-Bayesian approximation (see [mediation::mediate()]).
#' Default is 1000.
#' @param boot_ci_type 	a character string indicating the type of bootstrap
#' confidence intervals for when \code{boot = TRUE}. If \code{"bca"},
#' bias-corrected and accelerated (BCa) confidence intervals will be estimated.
#' If \code{"perc"}, percentile confidence intervals will be estimated
#' (see [mediation::mediate()]). Default is "bca".
#' @param ci_level the designated confidence level. Default 0.95.
#' @param fused a logical variable for whether the fused LASSO should be used
#' instead of the ordinary LASSO. Default is \code{FALSE}.
#' @param gamma a numeric variable \code{>=0} indicating the ratio of the
#' ordinary LASSO penalty to the fusion penalty (see [genlasso::genlasso()]).
#' Ignored if \code{fused = FALSE}. Default is 0.
#' @param eps numeric variable indicating the multiplier for the ridge penalty
#' in case \code{X} is rank deficient (see  [genlasso::genlasso()]).
#' Default is \code{1e-4}.
#' @param maxsteps an integer specifying the maximum number of steps for the
#' algorithm before termination (see [genlasso::genlasso()]). Default
#' is 2000.
#'
#' @return A list containing:
#' \itemize{
#'     \item{loadings: }{a matrix of the sparse PC loadings.}
#'     \item{pcs: }{a matrix of the PCs.}
#'     \item{var_explained: the cumulative proportion of variance explained by the PCs.}
#'     \item{contributions: }{a data frame containing the estimates, confidence
#'     intervals, and p-values of the mediation contributions.}
#'     \item{effects: }{a data frame containing the estimated direct, global
#'     mediation, and total effects}
#' }
#'
#' @import MASS
#' @import mediation
#' @import Matrix
#' @import igraph
#' @import genlasso
#' @import sandwich
#'
#' @references Zhao, Y., Lindquist, M. A. & Caffo, B. S. Sparse principal
#' component based high-dimensional mediation analysis. Comput. Stat.
#' Data Anal. 142, 106835 (2020).
#'
#' @export
#'

spcma <-
  function(A, M, Y, var_per = 0.8, adaptive = FALSE, n_pc = NULL, sims = 1000,
           boot_ci_type = "bca", ci_level = 0.95, fused = FALSE,  gamma = 0,
           eps = 1e-4, maxsteps = 2000){

    #Create column names if absent
    if (is.null(colnames(M))) colnames(M) <- paste0("m", 1:p)

    #Obtain sparse principal component loadings
    spca_out <-
      spca_loadings(A, M, adaptive = adaptive, fused = fused, var.per = var_per,
                    n.pc = n_pc, gamma = gamma, eps = eps, trace = FALSE,
                    maxsteps = maxsteps, lambda.tune = "R2",
                    per.jump = 0.70)

    #Organize loadings
    n_pc <- ncol(spca_out$U)
    loadings <- spca_out$W
    rownames(loadings) <- colnames(M)
    colnames(loadings) <- paste0("u", 1:n_pc)

    #Compute and decorrelate sparse PCs
    pcs_correlated <- M %*% loadings
    pcs <- decorrelate_MX(pcs_correlated, A)
    colnames(pcs) <- paste0("spc", 1:n_pc)

    #Run marginal mediation on sparse PCs
    spcma_out <-
      mediate_multiple(A, pcs, Y, sims = sims, boot.ci.type = boot_ci_type,
                       conf.level = ci_level)

    #Organize mediation contributions
    contributions <-
      with(
        spcma_out,
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
    percentiles <- round(c((1 - conf.level) /2, 1 - (1 - conf.level) / 2) * 100,1)
    ci_names <-  paste0("ab_cl_", percentiles,"%")
    colnames(contributions)[6:7] <- ci_names

    #Organize effects
    effects <- matrix(NA, 3, 6)
    effects[, 1] <- c("indirect","direct","total")
    effects[2, 2:5] <- spcma_out$DE[1, c(1, 2, 4:5, 3)]
    effects[1, 2:5] <- spcma_out$IE.total[1, c(1, 2, 4:5, 3)]
    effects[3, 2:5] <- spcma_out$TE[1, c(1, 2, 4:5, 3)]
    ci_names <- paste0("cl_",percentiles,"%")
    colnames(effects) <- c("effect","estimate","se",ci_names,"ab_pv")

    #Return output
    out <-
      list(
        loadings = as.data.frame(loadings),
        pcs = as.data.frame(pcs),
        var_explained = spca_out$var.spc,
        contributions = contributions,
        effects = as.data.frame(effects)
      )

    return(out)
  }

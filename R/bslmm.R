#' Bayesian Sparse Linear Mixed Model for Mediation Analysis with High-Dimensional
#' Mediators
#'
#' @description
#'
#' @param A length \code{n} numeric vector containing exposure variable
#' @param M \code{n x p} numeric matrix of high-dimensional mediators.
#' @param Y length \code{n} numeric vector containing continuous outcome variable.
#' @param C1 optional numeric matrix of covariates to include in the outcome model.
#' @param C2 optional numeric matrix of covariates to include in the mediator
#' model. Default is \code{C1}.
#' @param burnin number of MCMC draws prior to sampling.
#' @param ndraws number of MCMC draws after burn-in.
#' @param ci_level the desired credible interval level. Default 0.95.
#' @param weights optional numeric vector of observation weights.
#' @param k shape parameter for inverse gamma priors. Default 2.
#' @param lm0 scale parameter for inverse gamma prior on small normal components.
#' Default \code{1e-4}. If \code{k=2}, this is the prior mean on the variance of
#' the small normal components.
#' @param lm1 scale parameter for inverse gamma prior on large normal components.
#' Default 1. If \code{k=2}, this is the prior mean on the variance of the
#' large normal components.
#' @param l scale parameter for the other inverse gamma priors.
#'
#' @return A list containing:
#' \itemize{
#'     \item{contributions: }{a data frame containing the estimates, Bayesian
#'     credible intervals, and posterior inclusion probabilities (PIPs) of the
#'     mediation contributions.}
#'     \item{effects: }{a data frame containing the estimated direct, global
#'     mediation, and total effects}
#' }
#'
#' @import bama
#'
#' @references Song, Y. et al. Bayesian shrinkage estimation of high dimensional
#' causal mediation effects in omics studies. Biometrics 76, 700–710 (2020).
#'
#' @export
#'
#' @examples
bslmm <- function(A, M, Y, C1 = NULL, C2 = C1, burnin = 30000, ndraws = 5000,
                  ci_level = 0.95, weights = NULL, inits = NULL, k = 2,
                  lm0 = 1e-4, lm1 = 1, l = 1){

  p <- ncol(M)
  n <- nrow(M)

  if(ci_level >= 1 | ci_level <= 0){
    stop("Confidence level should be between 0 and 1.")
  }

  if(is.null(colnames(M))){
    colnames(M) <- paste0(m,1:p)
  }

  if(is.null(C1)){
    C1 <- matrix(1,n,1)
  }

  if(is.null(C2)){
    C2 <- matrix(1,n,1)
  }

  controls <- list(k = k, lm0 = lm0, lm1 = lm1, l = l)

  bama_out <- bama(Y, A, M, C1 = C1, C2 = C2, method = "BSLMM",
                   burnin = burnin, ndraws = ndraws + burnin,
                   weights = weights, control = controls)

  #Organize mediation contributions
  ab <- with(bama_out, alpha.a * beta.m)
  percentiles <- c((1 - ci_level) /2, 1 - (1 - ci_level) / 2)
  contributions <-
    with(
      bama_out,
      data.frame(
        mediator = colnames(M),
        alpha = colMeans(alpha.a),
        beta = colMeans(beta.m),
        alpha_beta = colMeans(ab),
        ab_posterior_sd = apply(ab, 2, sd),
        cl1 = apply(ab, 2, quantile, probs = percentiles[1]),
        cl2 = apply(ab, 2, quantile, probs = percentiles[2]),
        ab_pip = colMeans(alphaa_member * betam_member),

      )
    )

  ci_names <-  paste0("ab_cl_", round(percentiles * 100, 1), "%")
  colnames(contributions)[6:7] <- ci_names

  #Organize mediation effects
  effects <- matrix(NA, 3, 5)
  effects[,  1] <- c("indirect","direct","total")
  ci_names <-  paste0("cl_", round(percentiles * 100, 1), "%")
  colnames(effects) <- c("effect","estimate","posterior_sd",ci_names,"pvalue")

  gie <- rowSums(ab) #global indirect effects
  effects[1, 2] <- mean(gie)
  effects[1, 3] <- sd(gie)
  effects[1, 4:5] <- quantile(gie,percentiles)

  de <- bama_out$beta.a #direct effects
  effects[2, 2] <- mean(de)
  effects[2, 3] <- sd(de)
  effects[2, 4:5] <- quantile(de,percentiles)

  te <- gie + de #total effects
  effects[3, 2] <- mean(te)
  effects[3, 3] <- sd(te)
  effects[3, 4:5] <- quantile(te,percentiles)

  output <-
    list(
      contributions = contributions,
      effects = as.data.frame(effects)
    )

  return(output)
}





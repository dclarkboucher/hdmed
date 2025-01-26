#' Bayesian Sparse Linear Mixed Model
#'
#' @description
#' \code{mediate_bslmm} fits the Bayesian sparse linear mixed model proposed by
#' Song et al. (2020) for high-dimensional mediation analysis, estimating the
#' mediation contributions of potential mediators.
#'
#' @param A length \code{n} numeric vector containing exposure variable
#' @param M \code{n x p} numeric matrix of high-dimensional mediators.
#' @param Y length \code{n} numeric vector containing continuous outcome variable.
#' @param C1 optional numeric matrix of covariates to include in the outcome model. Do
#' not include an intercept term.
#' @param C2 optional numeric matrix of covariates to include in the mediator
#' model. Default is \code{C1}. Do not include an intercept term.
#' @param burnin number of MCMC draws prior to sampling.
#' @param ndraws number of MCMC draws after burn-in.
#' @param ci_level the desired credible interval level. Default is 0.95.
#' @param weights optional numeric vector of observation weights.
#' @param k shape parameter for the inverse gamma priors. Default is 2.
#' @param lm0 scale parameter for the inverse gamma prior on the variance of the
#' smaller-variance normal components. Default is \code{1e-4}. If \code{k=2},
#' this parameter equals the prior mean on the smaller normal variance.
#' @param lm1 scale parameter for the inverse gamma prior on the variance of the
#' larger-variance components of \code{beta_m}. Default is 1. If \code{k=2},
#' this parameter equals the prior mean on the larger normal variance of the
#' mediator-outcome associations.
#' @param lma1 scale parameter for the inverse gamma prior on the variance of the
#' larger-variance components of \code{alpha_a}. Default is 1. If \code{k=2},
#' this parameter equals the prior mean on the larger normal variance of the
#' exposure-mediator associations.
#' @param l scale parameter for the other inverse gamma priors.
#'
#' @details
#' \code{mediate_bslmm} is a wrapper function for the "BSLMM" option from [bama::bama()],
#' which fits a Bayesian sparse linear mixed model for performing mediation
#' analysis with high-dimensional mediators. The model assumes that
#' the mediator-outcome associations (\eqn{\beta_m}) and the exposure-mediator
#' associations (\eqn{\alpha_a}) independently follow a mixture of small-variance
#' and high-variance normal distributions, and that if a mediator \eqn{M_j} has both
#' \eqn{(\beta_m)_j} and \eqn{(\alpha_a)_j} belonging to the larger-variance distribution,
#' it has a notably large mediation contribution compared to the others. The
#' posterior inclusion probability (PIP) of belonging to both larger-variance
#' distributions is reported for each mediator as \code{ab_pip}.
#'
#' @return A list containing:
#'
#' * `contributions`: a data frame containing the estimates, Bayesian credible
#' intervals, and posterior inclusion probabilities of the mediation contributions
#'
#' * `effects`: a data frame containing the estimated direct, global mediation,
#' and total effects.
#'
#'
#' @importFrom bama bama
#'
#'
#' @references Song, Y. et al. Bayesian shrinkage estimation of high dimensional
#' causal mediation effects in omics studies. Biometrics 76, 700-710 (2020).
#'
#' @examples
#' A <- med_dat$A
#' M <- med_dat$M
#' Y <- med_dat$Y
#'
#' # Toy example with small burnin and ndraws
#' out <- mediate_bslmm(A, M, Y, burnin = 100, ndraws = 10)
#' out$effects
#' head(out$contributions)
#'
#'
#' @export
#'
mediate_bslmm <- function(A, M, Y, C1 = NULL, C2 = C1, burnin = 30000, ndraws = 5000,
                  ci_level = 0.95, weights = NULL, k = 2,
                  lm0 = 1e-4, lm1 = 1, lma1 = 1, l = 1){

  p <- ncol(M)
  n <- nrow(M)

  #Check A, M, Y
  if(is.data.frame(M)) M <- as.matrix(M)
  if(!is.numeric(A) | !is.vector(A)) stop("A must be numeric vector.")
  if(!is.numeric(M) | !is.matrix(M)) stop("M must be numeric matrix.")
  if(!is.numeric(Y) | !is.vector(Y)) stop("Y must be numeric vector.")
  if(is.null(colnames(M))){
    colnames(M) <- paste0("m",1:p)
  }

  #Check covariates
  if(!is.null(C1)){
    if(!is.numeric(C1) | !is.matrix(C1)){
      stop("C1 and C2 should be numeric matrices when specified.")
    }

    C1 <- as.data.frame(C1)
    int <- data.frame(int = rep(1,n))
    C1 <- cbind(int, C1)

  }

  if(!is.null(C2)){
    if(!is.numeric(C2) | !is.matrix(C2)){
      stop("C1 and C2 should be numeric matrices when specified.")
    }

    C2 <- as.data.frame(C2)
    int <- data.frame(int = rep(1,n))
    C2 <- cbind(int, C2)
  }

  if(ci_level >= 1 | ci_level <= 0){
    stop("Confidence level should be between 0 and 1.")
  }

  if(is.null(colnames(M))){
    colnames(M) <- paste0("m",1:p)
  }

  if(is.null(C1)){
    C1 <- matrix(1,n,1)
  }

  if(is.null(C2)){
    C2 <- matrix(1,n,1)
  }

  controls <- list(k = k, lm0 = lm0, lm1 = lm1, lma1 = lma1, l = l)

  bama_out <- bama::bama(Y, A, M, C1 = C1, C2 = C2, method = "BSLMM",
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
        ab_pip = colMeans(r1 * r3),
        row.names = NULL
      )
    )

  ci_names <-  paste0("ab_cl_", round(percentiles * 100, 1), "%")
  colnames(contributions)[6:7] <- ci_names

  #Organize mediation effects
  effects <- matrix(NA, 3, 5)
  effects[,  1] <- c("indirect","direct","total")
  ci_names <-  paste0("cl_", round(percentiles * 100, 1), "%")
  colnames(effects) <- c("effect","estimate","posterior_sd",ci_names)

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





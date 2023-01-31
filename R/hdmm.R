#' High-dimensional Multivariate Mediation Analysis with
#' Principal Directions of Mediation
#'
#' \code{mediate_hdmm} estimates the first  "direction of mediation" in the
#' causal mediation mechanism of an exposure \code{A}, an outcome \code{Y}, and high-dimensional
#' mediators \code{M} as proposed by Chén et al. (2018).
#'
#' @param A numeric vector containing exposure variable.
#' @param M numeric matrix of high-dimensional mediators. It is not recommended
#' to supply a matrix with more mediators than observations.
#' @param Y numeric vector containing continuous outcome variable.
#' @param sims number of Monte Carlo draws for nonparametric bootstrap or
#' quasi-Bayesian approximation. See \code{\link{mediate}}. Default is 1000.
#' @param boot_ci_type 	a character string indicating the type of bootstrap
#' confidence intervals for when \code{boot = TRUE}. If \code{"bca"},
#' bias-corrected and accelerated (BCa) confidence intervals will be estimated.
#' If \code{"perc"}, percentile confidence intervals will be estimated
#' (see [mediation::mediate()]). Default is "bca".
#' @param ci_level the designated confidence level. Default 0.95.
#' @param tol tolerance. Default 10^-5.
#' @param theta numeric vector of length 5 describing starting value of pathway
#' coefficients. Default is a vector of 1's.
#' @param w1 numeric vector of the same length of \code{A} specifying PDM
#' starting values. Default is a vector of 1's.
#' @param interval numeric vector proportional to the intervals from where the
#' smoothing parameter is searched. Default is 10^6.
#' @param step numerical number specifying step width for smoothing parameter
#' search. Default is 10^4.
#' @param imax integer specifying the maximum number of iterations allowed.
#' Default is 100.
#'
#' @details
#' HDMM provides latent variable approach to high-dimensional mediation analysis.
#' The function \code{mediate_hdmm} uses a likelihood-based approach to compute
#' principal directions of mediation (PDMs), which are loading weights used to linearly
#' combine the inputted mediators to form a single, latent variable that replaces
#' the original mediators in the analysis. Though HDMM cannot be used to estimate
#' the global mediation effect or the contributions of specific mediators, it can
#' still can be useful for inferring whether there is mediation occurring through
#' the set of mediators as a joint system. See the provided reference for more
#' details.
#'
#'
#' @return A list containing:
#' \itemize{
#'     \item{pdm: }{the first direction of mediation by which mediators are
#'     weighted.}
#'     \item{mediator: }{the latent mediator corresponding to the first direction
#'     of mediation.}
#'     \item{effects: }{a data frame containing the estimates, confidence
#'     intervals, and p-values of the mediation effects.}
#' }
#'
#'
#'
#' @references Chén, O. Y. et al. High-dimensional multivariate mediation with
#' application to neuroimaging data. Biostatistics (Oxford, England) vol. 19 (2018).
#'
#' @source <https://github.com/oliverychen/PDM>
#'
#'
#' @export
#'
#' @examples
#' A <- med_dat$A
#' M <- med_dat$M
#' Y <- med_dat$Y
#' out <- mediate_hdmm(A, M, Y, sims = 50)
#' out$effects
#'

mediate_hdmm <- function(A, M, Y, sims = 1000, boot_ci_type = "bca",
                 ci_level = 0.95, tol = 10^-5, theta = rep(1,5),
                 w1 = rep(1,ncol(M)), interval = 10^6, step = 10^4,
                 imax = 100){




  #Obtain first PDM and corresponding latent mediator
  pdm1 <- as.vector(PDM_1(A, Y, M, imax, tol, theta, w1, interval, step)$w1)
  m <- as.vector(M %*% as.matrix(pdm1, ncol = 1)) #latent mediator


  #Fit mediation models
  fit.total <- lm(Y ~ A)
  fit.m <- lm(m ~ A)
  fit.y <- lm(Y ~ A + m)

  #Mediation analysis
  med <-
    suppressMessages(
    mediation::mediate(fit.m, fit.y, treat = "A", mediator = "m", sims = sims, boot = T,
            boot.ci.type = boot_ci_type, conf.level = ci_level))

  #Create empty results table
  results <- matrix(NA, nrow = 3, ncol = 6)
  percentiles <- round(c((1 - ci_level) /2, 1 - (1 - ci_level) / 2) * 100,1)
  ci_names <-  paste0("cl_", percentiles,"%")
  colnames(results) <- c("effect","estimate","se",ci_names,"pv")

  #Fill in results
  results[1:3, 1] <- c("indirect","direct","total")
  results[1, 2:6] <- c(med$d1, sd(med$d1.sims), med$d1.ci, med$d1.p) #indirect effect
  results[2, 2:6] <- c(med$z1, sd(med$z1.sims), med$z1.ci, med$z1.p) #direct effect
  results[3, 2:6] <- c(med$tau.coef, sd(med$tau.sims), med$tau.ci, med$tau.p) #total effect

  out <- list(pdm = pdm1,
              mediator = m,
              effects = as.data.frame(results)
              )

  return(out)


}


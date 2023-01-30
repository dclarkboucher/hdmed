#' Mediation Analysis via Fixed Effects Model
#'
#' @description \code{mediate_medfix} fits a high-dimensional mediation model with
#' the adaptive LASSO approach as proposed by Zhang (2021) for the
#' special case of MedFix that there is only one exposure variable.
#'
#'
#' @param A length \code{n} numeric vector representing the exposure variable
#' @param M \code{n x p} numeric matrix of high-dimensional mediators.
#' @param Y length \code{n} numeric vector representing the continuous outcome variable.
#' @param C1 optional numeric matrix of covariates to include in the outcome model.
#' Default is \code{NULL}.
#' @param C2 optional numeric matrix of covariates to include in the mediator
#' model. Default is \code{C1}.
#' @param nlambda number of lambdas attempted in the adaptive LASSO.
#' See [gcdnet::cv.gcdnet()]. The specific sequence of \code{lambda}s is chosen
#' by the \code{cv.gcdnet()} function.
#' @param nlambda2 number of \code{lambda2}s attempted in the initial elastic net
#' used for computing adaptive weights prior to adaptive LASSO. Default is \code{50}.
#' See [gcdnet::cv.gcdnet()] for details on elastic net regression. If \code{0},
#' \code{lambda2=0} is fed to \code{cv.gcdnet()} and the initial fit is based on
#' standard LASSO, not elastic net. If not \code{0)}, the specific sequence of
#' \code{lambda2}s is given by \code{exp(seq(1e-4,0.02,length.out=nlambda2))}, and
#' the \code{lambda2} with the least cross-validated error is chosen.
#' @param nfolds number of folds for cross-validation. See [gcdnet::cv.gcdnet()].
#' Default is \code{10}.
#' @param seed numeric random seed.
#'
#'
#' @import gcdnet
#'
#' @details
#' MedFix performs mediation analysis when there are multiple mediators by
#' applying adaptive LASSO to the outcome model. In order to fit the adaptive LASSO,
#' we first obtain an initial model fit using either LASSO (which deploys
#' the L1 penalty) or elastic net (which deploys both the L1 and L2 penalties)
#' depending on the provided \code{nlambda2}. Estimates from this fit are then used
#' to compute the adaptive weights used in the adaptive LASSO. Once the final
#' adaptive LASSO estimates (\eqn{\beta_m}) are obtained for the outcome model, estimates for
#' the \code{p} mediator models (\eqn{\alpha_a}) are obtained by linear regression. The
#' mediation contributions are computed as \eqn{\alpha_a} times \eqn{\beta_m}, and the p-value
#' is taken as the maximum of the \eqn{\alpha_a} and \eqn{beta_m} p-values. Last, the
#' global indirect effect is estimated by summing the mediation contributions, and the direct
#' effect is estimated by subtracting the global indirect effect from an estimate of
#' the total effect. This function is specific to applying MedFix to the special
#' case that there is only one exposure; for details on how to apply
#' MedFix when the exposures are high-dimensional, as proposed by the
#' authors, see the supplemental files of the referenced manuscript.
#'
#'
#' @return A list containing:
#' \itemize{
#'     \item{contributions: }{a data frame containing the estimates and p-values
#'     of the mediation contributions}
#'     \item{effects: }{a data frame containing the estimated direct, global
#'     mediation, and total effects}
#' }
#'
#' @references
#' Zhang, Q. High-Dimensional Mediation Analysis with Applications to Causal
#' Gene Identification. Statistics in biosciences (2021)
#'
#' @examples
#' data("med_dat")
#' A <- med_dat$Y
#' M <- med_dat$M
#' Y <- med_dat$Y
#'
#' out <- mediate_medfix(A, M, Y, nlambda = 10, nlambda2 = 5, seed = 1)
#'
#' @export
#'

mediate_medfix <- function(A, M, Y, C1 = NULL, C2 = C1,
                           nlambda = 100, nlambda2 = 50, nfolds = 10,
                           seed = 1){

  n <- nrow(M)
  p <- ncol(M)
  set.seed(seed)

  #Check A, M, Y
  if(is.data.frame(M)) M <- as.matrix(M)
  if(!is.numeric(A) | !is.vector(A)) stop("A must be numeric vector.")
  if(!is.numeric(M) | !is.matrix(M)) stop("M must be numeric matrix.")
  if(!is.numeric(Y) | !is.vector(Y)) stop("Y must be numeric vector.")
  if(is.null(colnames(M))){
    colnames(M) <- paste0("m",1:ncol(M))
  }

  #Check covariates
  if(!is.null(C1)){
    if(!is.numeric(C1) | !is.matrix(C1)){
      stop("C1 and C2 should be numeric matrices when specified.")
    }
  }

  if(!is.null(C2)){
    if(!is.numeric(C2) | !is.matrix(C2)){
      stop("C1 and C2 should be numeric matrices when specified.")
    }
  }

  #Implement MedFix
  med_out <- medfix(A, M, Y, C1, C2, nlambda = nlambda, nlambda2 = nlambda2)

  #Organize mediation contributions
  contributions <- med_out[, -6] #remove total effect

  #Organize effects
  te <- med_out$te[1]
  gie <- sum(med_out$alpha_beta)
  de <- te - gie
  effects <-
    data.frame(
      effect = c("indirect", "direct", "total"),
      estimate = c(gie,de,te)
    )

  output <-
    list(
      contributions = contributions,
      effects = effects
    )

  return(output)


}

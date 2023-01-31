#' High-Dimensional Mediation Analysis
#'
#' @description \code{mediate_hdma} fits a high-dimensional mediation model with
#' the de-biased LASSO approach as proposed by Gao et al. (2022),
#' estimating the mediation contributions of potential mediators.
#'
#' @param A length \code{n} numeric vector containing exposure variable
#' @param M \code{n x p} numeric matrix of high-dimensional mediators.
#' @param Y length \code{n} numeric vector containing continuous or binary outcome variable.
#' @param C1 optional numeric matrix of covariates to include in the outcome model.
#' @param C2 optional numeric matrix of covariates to include in the mediator model.
#' @param binary_y logical flag for whether \code{Y} should be interpreted as a
#' binary variable with 1/0 coding rather than as continuous. Default is \code{FALSE}.
#' @param n_include integer specifying the number of top markers from sure
#' independent screening to be included. Default is \code{NULL}, in which case
#' \code{n_include} will be either \code{ceiling(n/log(n))} if
#' \code{binary_Y = F}, or \code{ceiling(n/(2*log(n)))} if \code{binary_Y = T}.
#' If \code{n_include >= p}, all mediators are included with no screening. Note
#' that if \code{binary_y = F}, screening is performed based on the single-mediator
#' outcome model p-values, and if \code{binary_y = F}, screening is based on the
#' the mediator model p-values.
#' @param ... other arguments passed to [hdi::hdi()].
#'
#'
#' @details The first step in HDMA is to perform sure independence
#' screening (SIS) to choose the \code{n_include} mediators that are most
#' associated with the outcome (when Y is continuous) or the exposure
#' (when Y is binary), based on p-values from linear regression. The second step
#' is to fit the outcome model for the remaining mediators using de-sparsified
#' (A.K.A de-biased) LASSO, which as asymptotic properties allowing for
#' computation of p-values by the \code{hdi} package. HDMA then fits the
#' mediator models using linear regression among those mediators that have both
#' survived SIS (in step 1) and been identified by the LASSO (in step 2), obtaining
#' p-values for the mediation contributions by taking the maximum of the \eqn{\alpha_a}
#' and \eqn{\beta_m} p-values. The global indirect effect is estimated by summing the
#' mediation contributions, and the direct effect is estimated by subtracting
#' the global indirect effect from an estimate of the total effect. See References for
#' more detail.
#'
#' @return A list containing:
#' \itemize{
#'     \item{contributions: }{a data frame containing the estimates and p-values
#'     of the mediation contributions}
#'     \item{effects: }{a data frame containing the estimated direct, global
#'     mediation, and total effects}
#' }
#'
#'
#'
#' @references
#' Gao, Y. et al. Testing Mediation Effects in High-Dimensional
#' Epigenetic Studies. Front. Genet. 10, 1195 (2019).
#'
#' Fan, J. & Lv, J. Sure independence screening for ultrahigh dimensional
#' feature space. J. R. Stat. Soc. 70, 849–911 (2008).
#'
#' @source <https://github.com/YuzhaoGao/High-dimensional-mediation-analysis-R>
#'
#' @examples
#' A <- med_dat$A
#' M <- med_dat$M
#' Y <- med_dat$Y
#'
#' # Fit hdma with continuous outcomes
#' out <- mediate_hdma(A, M, Y)
#' head(out$contributions)
#' out$effects
#'
#' @export
#'

mediate_hdma <- function(A, M, Y, C1 = NULL, C2 = NULL, binary_y = F,
                         n_include = NULL,
                         ...){

  n <- nrow(M)
  p <- ncol(M)

  if (!is.null(n_include)){
    n_include <- min(n_include, p)
  }

  if (binary_y){

    family <- "binomial"
    if(!all(Y %in% c(0,1))){
      stop("Please code Y as 0 or 1 when invoking the binary_y flag.")
    }

  } else{
    family <- "gaussian"

  }

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
  }

  if(!is.null(C2)){
    if(!is.numeric(C2) | !is.matrix(C2)){
      stop("C1 and C2 should be numeric matrices when specified.")
    }

    C2 <- as.data.frame(C2)
  }

  #Fit HDMA
  hdma_out <- hdma(A, Y, M, COV.XM = C2, COV.MY = C1, family = family,
                   topN = n_include, method = "lasso", ...)

  if (nrow(hdma_out) == 0) return(NULL)

  #Organize mediation contributions
  contributions <- hdma_out[, 1:7]

  #Organize effects
  te <- hdma_out$te[1]
  gie <- sum(hdma_out$alpha_beta)
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




#' High-Dimensional Mediation Analysis
#'
#' @description \code{mediate_hima} fits a high-dimensional mediation model with
#' the minimax concave penalty as proposed by Zhang et al. (2016),
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
#' @param ... other arguments passed to \code{\link{hdi}}.
#'
#'
#' @details The first step in HIMA is to perform sure independence
#' screening (SIS) to choose the \code{n_include} mediators that are most
#' associated with the outcome (when Y is continuous) or the exposure
#' (when Y is binary), based on p-values from linear regression. The second step
#' is to fit the outcome model for the remaining mediators with the minimax
#' concave penalty. HIMA then fits the mediator models using linear regression
#' among those mediators that have both survived SIS (in step 1) and been
#' selected by the MCP (in step 2), which enables estimation of the mediation
#' contributions. The global indirect effect is estimated by summing these
#' contributions, and the direct effect is estimated by subtracting the global
#' indirect effect from an estimate of the total effect. We compute p-values for
#' the mediation contributions by taking the maximum of the \eqn{\alpha_a} and \eqn{\beta_m}
#' p-values, where the beta p-values are obtained via a second, unpenalized generalized
#' linear model containing only the mediators selected by the MCP. We include this
#' p-value computation so that our function replicates the behavior of [HIMA::hima()],
#' the function on which ours is based, but we caution that the beta p-values
#' may be over-optimistic due to double-dipping, since the mediators tested in
#' the unpenalized model are only those chosen by the penalized model. Note also
#' that the HIMA authors apply Bonferroni correction to the final, maxed p-values
#' to account for multiple testing, which we choose to leave up to the user. For
#' more information, see the R package \code{\link{HIMA}} along with the provided
#' reference.
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
#' @import glmnet
#' @import ncvreg
#' @import matrix
#' @import stats
#' @import foreach
#' @import iterators
#'
#'
#' @references
#' Zhang, H. et al. Estimating and testing high-dimensional mediation effects
#' in epigenetic studies. Bioinformatics 32, 3150–3154 (2016).
#'
#' @examples
#' data("med_dat")
#' A <- med_dat$A
#' M <- med_dat$M
#' Y <- med_dat$Y
#'
#' # Fit hima with continuous outcome
#' out <- mediate_hima(A, M, Y)
#' head(out$contributions)
#' out$effects
#'
#' # Fit hima with binary outcome
#' Y1 <- as.numeric(Y > mean(Y))
#' out1 <- mediate_hima(A, M, Y1, binary_y = T)
#' head(out1$contributions)
#' out1$effects
#'
#' @export
#'
#'
#'
mediate_hima <- function(A, M, Y, C1 = NULL, C2 = NULL, binary_y = F,
                         n_include = NULL,
                         ...){

  n <- nrow(M)
  p <- ncol(M)

  if (!is.null(n_include)){
    n_include <- min(n_include, p)
  }

  if (binary_y){

    family <- "binomial"
    if(!all(y %in% c(0,1))){
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

  #Fit HIMA
  hima_out <- hima(A, Y, M, COV.XM = C2, COV.MY = C1, Y.family = family,
                   M.family = "gaussian", penalty = "MCP", ...)

  if(nrow(hima_out) == 0) return(NULL)

  #Organize mediation contributions
  contributions <- hima_out[,1:7]

  #Organize effects
  te <- hima_out$te[1]
  gie <- sum(hima_out$alpha_beta)
  de <- te - gie
  effects <-
    data.frame(
      effect = c("indirect", "direct", "total"),
      estimate = c(gie, de, te)
    )

  output <-
    list(
      contributions = contributions,
      effects = effects
    )

  return(output)

}




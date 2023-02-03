#' High-Dimensional Linear Mediation Analysis
#'
#'
#' @description \code{mediate_hilma} applies high-dimensional linear mediation
#' analysis (HILMA) as proposed by Zhou et al. (2020).
#'
#' @param A length \code{n} numeric vector containing a single exposure variable
#' or size \code{n x q} numeric matrix containing multiple exposures. If a matrix,
#' should be low-dimensional (\code{q << n}).
#' @param M \code{n x p} numeric matrix of high-dimensional mediators.
#' @param Y length \code{n} numeric vector containing continuous outcome variable.
#' @param aic_tuning logical flag for whether to select the tuning parameter using
#' AIC. Default is \code{FALSE} as this was not the preferred approach by the
#' HILMA authors (see references for more detail).
#' @param nlambda number of candidate lambdas for AIC tuning. Default is 5. If
#' \code{aic_tuning=F} this parameter is ignored.
#' @param lambda_minmax_ratio ratio of the minimum lambda attempted in
#' AIC tuning to the maximum. If \code{aic_tuning=F}, ignored.
#' @param center logical flag for whether the variables should be centered. Default
#' is \code{TRUE}.
#'
#' @details
#' \code{mediate_hilma} is a wrapper function for the [freebird::hilma()] function,
#' which fits the "high-dimensional linear mediation analysis" model proposed by
#' Zhou et al. (2020) for mediation settings when there are high-dimensional
#' mediators and one or several exposures. The function returns estimates of
#' the direct effect, total effect, and global mediation effect, the last of which
#' is tested for statistical significance with a reported p-value. For additional
#' detail, see the attached reference as well as the [freebird::hilma()]
#' documentation.
#'
#'
#' @return A list containing, for each exposure, a data frame of
#' the estimated direct, total, and global mediation effects. A p-value is
#' provided for the global mediation effect.
#'
#' @importFrom freebird hilma
#'
#' @references Zhou, R. R., Wang, L. & Zhao, S. D. Estimation and inference for
#' the indirect effect in high-dimensional linear mediation models. Biometrika
#' 107, 573–589 (2020).
#'
#'
#' @examples
#' A <- med_dat$A
#' M <- med_dat$M
#' Y <- med_dat$Y
#'
#' # Implement HILMA with one exposure
#' out <- mediate_hilma(A, M, Y)
#' out$a1
#'
#'
#'
#' @export
#'
mediate_hilma <- function(A, M, Y, aic_tuning = F, nlambda = 5,
                          lambda_minmax_ratio = 0.1, center = T){

  n <- nrow(M)
  p <- ncol(M)

  #Check A, M, Y
  if(is.data.frame(M)) M <- as.matrix(M)
  if(!is.numeric(A)) stop("A must be numeric vector or matrix.")
  if(!is.numeric(M) | !is.matrix(M)) stop("M must be numeric matrix.")
  if(!is.numeric(Y) | !is.vector(Y)) stop("Y must be numeric vector.")
  if(is.null(colnames(M))){
    colnames(M) <- paste0("m",1:p)
  }

  if (is.vector(A)){
    A <- matrix(A, n, 1)
    colnames(A) <- "a1"
  }

  q <- ncol(A)
  if(is.null(colnames(A))){

    colnames(A) <- paste0("a",1:q)
  }

  nlambda <- as.integer(nlambda)
  if(aic_tuning & (is.na(nlambda) | nlambda <=0)){
    stop("nlambda must be a positive integer.")
  }

  if(lambda_minmax_ratio <=0 | lambda_minmax_ratio >= 1){
    stop("lambda_minmax_ratio should be between 0 and 1.")
  }

  if (!aic_tuning){

    hilma_out <- hilma(Y, M, A, mediation_setting = "incomplete",
                       tuning_method = "uniform", center = center,
                       n.lambda = 1)

  } else{

    hilma_out <- freebird::hilma(Y, M, A, mediation_setting = "incomplete",
                                 tuning_method = "aic", center = center,
                                 min.ratio = lambda_minmax_ratio, n.lambda = nlambda)


  }


  output <- lapply(1:q, extract_hilma, hilma_out = hilma_out)

  names(output) <- colnames(A)

  return(output)


}

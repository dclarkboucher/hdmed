#' High-Dimensional Linear Mediation Analysis
#'
#'
#' @description \code{mediate_hilma} applies high-dimensional linear mediation
#' analysis (HILMA) as proposed by Zhou et al. (2020).
#'
#' @param A length \code{n} numeric vector containing exposure variable or size
#' \code{n x q} numeric matrix containing multiple exposure variables.
#' @param M \code{n x p} numeric matrix of high-dimensional mediators.
#' @param Y length \code{n} numeric vector containing continuous outcome variable.
#' @param aic_tuning logical flag for whether to select the tuning parameter using
#' AIC. Default is \code{FALSE} as this was not the preferred approach by the
#' HILMA authors (see references for more detail).
#' @param nlambda number of candidate lambdas for AIC tuning. Default 5. If
#' \code{aic_tuning=F}, ignored.
#' @param lambda_minmax_ratio ratio of the minimum lambda attempted in
#' AIC tuning to the maximum. If \code{aic_tuning=F}, ignored.
#' @param center logical flag for whether the variables should be centered. Default
#' is \code{TRUE}.
#'
#' @return A list containing, for each exposure variable, a data frame of
#' the estimate direct, total, and global mediation effects. A p-value is
#' provided for the global mediation effect.
#'
#'
#' @importFrom freebird hilma
#'
#' @references Zhou, R. R., Wang, L. & Zhao, S. D. Estimation and inference for
#' the indirect effect in high-dimensional linear mediation models. Biometrika
#' 107, 573–589 (2020).
#'
#' @export
#'
#' @examples
mediate_hilma <- function(A, M, Y, aic_tuning = F, nlambda = 5,
                          lambda_minmax_ratio = 0.1, center = T){

  n <- nrow(M)
  p <- ncol(M)

  if (is.vector(A)){
    A <- matrix(A,n,1)
  }

  q <- ncol(A)

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

    hilma_out <- hilma(Y, M, A, mediation_setting = "incomplete",
                       tuning_method = "aic", center = center,
                       min.ratio = lambda_minmax_ratio, n.lambda = nlambda)


  }

  message(paste("lambda used:",hilma_out$lambda_used))

  output <- lapply(1:q, extract_hilma, hilma_out = hilma_out)

  if (!is.null(colnames(A))){
    names(output) <- colnames(A)

  } else{
    names(output) <- paste0("a",1:q)

  }

  return(output)


}

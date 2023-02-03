#' Latent Variable Mediation Analysis
#'
#' @description \code{mediate_lvma} fits a high-dimensional mediation model described
#' by Derkach et al. (2019), in which a small number of latent, unmeasured
#' mediators replace the original mediators in the model.
#'
#'
#' @param A length \code{n} numeric vector representing the exposure variable
#' @param M \code{n x p} numeric matrix of high-dimensional mediators.
#' @param Y length \code{n} numeric vector representing the continuous outcome variable.
#' @param q number of latent mediators
#' @param scale logical flag for whether the inputted mediators should be standardized
#' prior to the analysis. Default is \code{FALSE}, but \code{TRUE} may be worth
#' attempting in case of errors.
#' @param rhoLM numeric vector of candidate penalty parameters for the latent
#' mediator-mediator associations in the joint likelihood. Default is a short
#' toy sequence.
#' @param rhoEL numeric vector of candidate penalty parameters for the exposure-latent mediator
#' associations in the joint likelihood. Default is a short toy sequence.
#' @param rhoLY numeric vector of candidate penalty parameters for the
#' latent mediator-outcome associations in the joint likelihood. Default is a
#' short toy sequence.
#' @param imax integer specifying the maximum number of iterations allowed. Default is 5000.
#'
#' @details
#' LVMA is a latent variable mediation model which assumes, contrary to standard
#' assumptions, that the inputted set of candidate mediators do not
#' affect the outcome through the exposure on their own, but rather, occur as result
#' of latent, unmeasured mediators which themselves transmit effects from the
#' exposure to outcome. The required parameters for fitting this model are
#' \code{rhoLE}, a regularization parameter for effects of the latent mediators
#' on the inputted mediators; \code{rhoEL}, a regularization parameter for the
#' effects of the exposure on the latent mediators; and \code{rhoLY}, a
#' regularization parameter for the effects of the latent mediators on the
#' exposure. These parameters should ideally be supplied by the user as
#' vectors, so that each combination of the three parameters can be attempted
#' in the estimation. However, this can be intensely computation costly, and
#' for simplicity our default values are vectors of length 4, corresponding
#' to a 64 by 64 parameter grid. In practice, Derkach et al. use a much larger
#' grid with 5 values of \code{rhoLM} (ranging from 6 to 8.5), 40 values of
#' \code{rhoEY} (ranging from 0 to 40), and 40 values of \code{rhoLY} (ranging
#' from 0 to 75). Supplying longer parameter vectors makes the fit more flexible,
#' but more computationally costly, and to reliably implement LVMA on real data
#' one should use a larger parameter grid with parallel computation on a remote
#' computing cluster, as did the authors. For more information on the likelihood,
#' parameters, and mediation model, see the referenced article and/or its
#' supplement files.
#'
#' @return A list containing the selected models based on AIC, BIC, and EBIC
#' (recommended) as three sub-lists. The sub-lists include objects indicating
#' the penalty set that was used (\code{penalty}), the values of the chosen parameters
#' (e.g., \code{EBIC}), the exposure-latent mediator effects (\code{AL_effects}),
#' the latent mediator-mediator effects (\code{LM_effects}, a data frame), the
#' direct effect of the exposure on the outcome (\code{AY_direct_effect}), the
#' the latent mediator-outcome effects (\code{LY_effects}), and binary vector
#' indicating whether each mediator was determined to be active. Here, active
#' mediators are those which are associated with a latent mediator that itself
#' is associated with both A and Y.
#'
#'
#' @references
#' Derkach, A., Pfeiffer, R. M., Chen, T.-H. & Sampson, J. N. High dimensional
#' mediation analysis with latent variables. Biometrics 75, 745–756 (2019).
#'
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/30859548/}
#'
#' @examples
#' A <- med_dat$A
#' M <- med_dat$M
#' Y <- med_dat$Y
#'
#' # Perform latent variable mediation analsis with 4 latent mediators and print
#' # whether the original 20 mediators are "actively" related to mediation
#' out <- mediate_lvma(A, M, Y, q = 4, rhoLM = 2, rhoEL = 2, rhoLY = 2, imax = 50)
#' table(out$EBIC_out$mediator_active)
#'
#' @export
#'

mediate_lvma <- function(A, M, Y, q, rhoLM, rhoEL, rhoLY, scale = T, imax = 5000){

  n <- nrow(M)
  p <- ncol(M)

  #Check A, M, Y
  if(is.data.frame(M)) M <- as.matrix(M)
  if(!is.numeric(A) | !is.vector(A)) stop("A must be numeric vector.")
  if(!is.numeric(M) | !is.matrix(M)) stop("M must be numeric matrix.")
  if(!is.numeric(Y) | !is.vector(Y)) stop("Y must be numeric vector.")

  #LVMA
  lvma_out <- lvma(
    A,
    M,
    Y,
    q = q,
    rhoL = rhoLM,
    rhoE = rhoEL,
    rhoY = rhoLY,
    imax = imax
  )

  results <- c()
  for(metric in c("EBIC", "BIC", "AIC")){
    results[[paste0(metric,"_out")]] <-
      process_lvma(
        M = M,
        out = lvma_out,
        type = metric,
        rhoLM = rhoLM,
        rhoEL = rhoEL,
        rhoLY = rhoLY
      )
  }

  return(results)


}

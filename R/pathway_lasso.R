#' Pathway LASSO for Mediation Analysis with High-Dimensional Mediators
#'
#' @description \code{mediate_plasso} fits a high-dimensional mediation model
#' with the penalized likelihood described by Zhao and Luo (2022), estimating the
#' mediation contributions of each site.
#'
#' @param A numeric vector containing exposure variable.
#' @param M numeric matrix of high-dimensional mediators.
#' @param Y numeric vector containing continuous outcome variable.
#' @param lambdas numeric vector of tuning parameters \code{lambda} for
#' which model is fitted. Default is a vector of length 50 ranging from 10^-5
#' to 10^4, with more density in the lower range.
#' @param select_lambda logical flag indicating whether to conduct a tuning
#' parameter selection using the Variable Selection Stability Criterion
#' described by Sun et al. (2013). Default is \code{FALSE}.
#' @param vss_rep if \code{select_lambda} is \code{TRUE}, number of VSSC
#' replications.
#' @param vss_cutoff if \code{select_lambda} is \code{TRUE}, cutoff used for
#' VSSC. Default is 0.1.
#' @param omega_ratio ratio of the \code{omega} parameter to \code{lambda}
#' parameter in the likelihood penalty. Default is 1.
#' @param phi value of the \code{phi} parameter in the likelihood penalty.
#' Default is 2. Cannot be less than 1/2.
#' @param maxit the maximum number of iterations. Default is 5000.
#' @param tol convergence tolerance. Default is 10^-6.
#'
#'
#' @details
#' Pathway LASSO fits a high-dimensional mediation model with a likelihood
#' that directly penalizes the mediator-outcome effects, exposure-mediator
#' effects, and mediation contributions (i.e., the mediation "pathways").
#' The shrinkage of the model is determined by three parameters--\code{phi},
#' \code{omega}, and \code{lambda}--with higher values tending to result in more
#' sparsity. Mediation results are returned for every unique value in the
#' inputted \code{lambdas} argument, in increasing order, excluding potentially
#' any for which the estimation was unsuccessful. For details on the exact
#' likelihood, see the first reference.
#'
#' When implementing multiple lambdas (i.e., using either the default argument
#' to \code{lambdas} or specifying one's own vector), there is an option to
#' perform a tuning parameter selection using the Variable Selection Stability
#' Criterion proposed by Sun et al. (2013), which chooses the parameter for
#' which the variable selection is most stable in repeated data-splitting.
#' However, implementing this can be computationally costly, especially with a
#' long list of tuning parameters, since it involves re-fitting pathway LASSO
#' many times.
#'
#' @return A list containing:
#' \itemize{
#'     \item{lambdas: }{the \code{lambda}s attempted in the same order as the objects in \code{all_fits}.}
#'     \item{all_fits: }{a list containing, for each \code{lambda}, a \code{data.frame} of the estimated mediation effects.}
#'     \item{chosen_lambda: }{if \code{select_lambda} is \code{TRUE}, the \code{lambda} chosen by VSSC.}
#'     \item{chosen_fit: }{if \code{select_lambda} is \code{TRUE}, the fit corresponding to the chosen \code{lambda}.}
#'     \item{vss: }{if \code{select_lambda} is\code{TRUE}, a \code{data.frame} containing the variable selection stabilities.}
#' }
#'
#'
#' @import MASS
#'
#'
#' @references Zhao, Y. & Luo, X. Pathway LASSO: pathway estimation and
#' selection with high-dimensional mediators. Stat. Interface 15, 39-50 (2022).
#'
#' Sun, W., Wang, J. & Fang, Y. Consistent selection of tuning parameters via
#' variable selection stability. J. Mach. Learn. Res. 14, 3419-3440 (2013).
#'
#' @source \url{https://github.com/zhaoyi1026/PathwayLasso}
#'
#' @examples
#' A <- med_dat$A
#' M <- med_dat$M[,1:8]
#' Y <- med_dat$Y
#' # fit pathway LASSO for two tuning parameters and retrieve their fits
#' out <- mediate_plasso(A, M, Y, lambdas = c(10^-3, 10^-2), tol = 1e-4)
#' head(out$all_fits$lambda1)
#' head(out$all_fits$lambda2)
#'
#'
#' @export


mediate_plasso <- function(A, M, Y, lambdas = NULL, select_lambda = FALSE,
                           vss_rep = 5,vss_cutoff = 0.1, omega_ratio = 1,
                           phi = 2, maxit = 5000, tol = 1e-6){

  p <- ncol(M)

  # Check A, M, Y
  if (is.data.frame(M)) M <- as.matrix(M)
  if (!is.numeric(A) | !is.vector(A)) stop("A must be numeric vector.")
  if (!is.numeric(M) | !is.matrix(M)) stop("M must be numeric matrix.")
  if (!is.numeric(Y) | !is.vector(Y)) stop("Y must be numeric vector.")
  if (is.null(colnames(M))){
    colnames(M) <- paste0("m",1:p)
  }else if(any(c("R","M") %in% colnames(M))){ #Avoid colname issues...
    warning("Mediator names overwritten. Avoid naming mediators 'R' or 'Z'")
    colnames(M) <- paste0("m",1:k)
  }

  if(omega_ratio < 0){
    stop("Omega ratio should be at least 0.")
  }

  if(phi < 0.5){
    stop("Phi must be at least 1/2 for the optimization problem to be convex.")
  }

  Z <- A
  R <- Y

  k <- ncol(M)
  dd0 <- data.frame(Z = Z, M, R = R)
  Sigma10<-diag(rep(1, k))
  Sigma20<-matrix(1, 1, 1)

  # Record SDs and standardize
  sd.Z <- sd(Z)
  sd.M <- apply(M, 2, sd)
  sd.R <- sd(R)

  Z <- scale(Z)
  M <- scale(M)
  R <- scale(R)
  dd <- data.frame(Z=Z, M, R=R)

  # Lambda
  if(is.null(lambdas)){
    lambdas <- c(10^c(seq(-5,-3,length.out=10),
                      seq(-3,0,length.out=26)[-1],
                      seq(0,2,length.out=11)[-1],
                      seq(2,4,length.out=6)[-1])) # lambda values

    nlambda <- length(lambdas)

  }else{
    lambdas <- unique(lambdas)
    nlambda <- length(lambdas)
    if(is.unsorted(lambdas)){
      lambdas <- sort(lambdas)
      message("Sorting lambdas increasingly")
    }
  }

  names(lambdas) <- paste0("lambda", 1:nlambda)

  # A few other parameters
  rho <- 1
  thred <- 1e-6
  thred2 <- 1e-3
  zero.cutoff <- 1e-3

  # Empty result objects
  AB.est  <- matrix(NA, k, length(lambdas))
  A.est <- AB.est
  B.est <- AB.est
  C.est <- rep(NA, length(lambdas))
  fit_succeeded <- c()

  message("Fitting pathway LASSO...")
  for (i in nlambda:1){

    out <- NULL

    # If we have no fits yet, fit the model naively.
    # Otherwise, use the last successful fit as a burn-in
    if(!any(fit_succeeded)){
      try(out<-mediation_net_ADMM_NC(Z,M,R,lambda=lambdas[i],
                                     omega=omega_ratio*lambdas[i],
                                     phi=phi,Phi1=NULL,Phi2=NULL,
                                     rho=rho,rho.increase=FALSE,
                                     tol=tol,max.itr=maxit,thred=thred,
                                     Sigma1=Sigma10,Sigma2=Sigma20,
                                     trace=FALSE))

    }else{

      last_fit <- rev(which(fit_succeeded))[1]
      which_burnin <- nlambda + 1 - last_fit # which prior fit to use as burn-in

      try(out<-mediation_net_ADMM_NC(Z,M,R,lambda=lambdas[i],
                                     omega=omega_ratio*lambdas[i],
                                     phi=phi,Phi1=NULL,Phi2=NULL,
                                     rho=rho,rho.increase=FALSE,
                                     tol=tol,max.itr=maxit,
                                     thred=thred,Sigma1=Sigma10,
                                     Sigma2=Sigma20,trace=FALSE,
                                     Theta0=matrix(c(1,A.est[,which_burnin]*(sd.Z/sd.M)),nrow=1),
                                     D0=matrix(c(C.est[which_burnin]*(sd.Z/sd.R),B.est[,which_burnin]*(sd.M/sd.R)),ncol=1),
                                     alpha0=matrix(c(1,A.est[,which_burnin]*(sd.Z/sd.M)),nrow=1),
                                     beta0=matrix(c(C.est[which_burnin]*(sd.Z/sd.R),B.est[,which_burnin]*(sd.M/sd.R)),ncol=1)))


    }

    # Did the algorithm succeed?
    succeeded <- !is.null(out)
    fit_succeeded <- append(fit_succeeded, succeeded)

    if(succeeded){

      B.est[, i] <- out$B * (sd.R / sd.M)
      C.est[i] <- out$C * (sd.R / sd.Z)
      A.est[, i] <- out$A * (sd.M / sd.Z)
      AB.est[, i] <- A.est[, i] * B.est[, i]

    }

  }

  # Organize output
  fit_succeeded <- rev(fit_succeeded)
  if(!any(fit_succeeded)){
    message("Pathway LASSO failed for all attempted lambdas.")
    return(NULL)
  }

  result_list <- c()
  lambdas1 <- lambdas[fit_succeeded] # Succeeded lambdas
  A.est <- as.matrix(A.est[, fit_succeeded])
  B.est <- as.matrix(B.est[, fit_succeeded])
  AB.est <- A.est * B.est
  C.est <- C.est[fit_succeeded]

  nlambda1 <- length(lambdas1)
  for(i in 1:nlambda1){
    result_list[[i]] <-
      data.frame(
        mediator = colnames(M),
        lambda = lambdas1[i],
        alpha = A.est[, i],
        beta = B.est[, i],
        alpha_beta =  AB.est[, i],
        direct_effect = C.est[i],
        global_indirect_effect = sum(AB.est[, i]),
        total_effect = sum(AB.est[, i]) + C.est[i],
        row.names = NULL
      )
  }

  names(result_list) <- names(lambdas1)

  output <-
    list(
      all_fits = result_list,
      lambdas = lambdas1
    )

  if(!select_lambda){
    return(output)
  }

  if(nlambda == 1){
    message("Only one lambda provided. Parameter selection ignored.")
    return(output)
  }

  message("All fits complete. Beginning tuning parameter selection by VSSC")
  vss_results <-
    mediation_net_ADMM_NC_KSC(Z,M,R,zero.cutoff=zero.cutoff,n.rep=vss_rep,
                              vss.cut=vss_cutoff,lambda=lambdas1,
                              omega=omega_ratio*lambdas1,
                              phi=phi,Phi1=NULL,Phi2=NULL,rho=rho,
                              rho.increase=FALSE,tol=tol,max.itr=maxit,
                              thred=thred, Sigma1=Sigma10,Sigma2=Sigma20,trace=FALSE,
                              Theta0=NULL,D0=NULL,alpha0=NULL,beta0=NULL)

  which_suggested_lambdas <- vss_results$lambda.idx
  chosen_lambda <- vss_results$lambda.est
  which_chosen <- which(lambdas1 == chosen_lambda)

  vss_out <-
    data.frame(
      index = names(lambdas1),
      lambda = lambdas1,
      mean_vss = vss_results$vss,
      meets_cutoff = (1:nlambda) %in% which_suggested_lambdas,
      chosen = lambdas1 == chosen_lambda, #minimum lambda meeting threshold is chosen
      row.names = NULL
    )

  output$chosen_lambda = lambdas1[which_chosen]
  output$chosen_fit = result_list[[which_chosen]]
  output$vss = vss_out

  return(output)

}



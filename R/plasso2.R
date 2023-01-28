
mediate_plasso <- function(A, M, Y, lambdas = NULL, select_lambda = F,
                           vss_rep = 5,vss_cutoff = 0.1, omega_ratio = 1,
                           phi = 2, maxit = 5000, tol = 1e-6){

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

  Z <- as.numeric(scale(Z))
  M <- as.matrix(scale(M))
  R <- as.numeric(scale(R))
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
  AB.est  <- matrix(NA, k, length(lambda))
  A.est <- AB.est
  B.est <- AB.est
  C.est <- rep(NA, length(lambda))
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

    }else{
      B.est[, i] <- NA
      C.est[i] <- NA
      A.est[, i] <- NA
      AB.est[, i] <- NA
    }

  }

  # Organize output
  fit_succeeded <- rev(fit_succeeded)
  if(!any(fit_succeeded)){
    message("Pathway LASSO failed for all attempted lambdas.")
    return(NULL)
  }

  result_list <- ls()
  lambdas1 <- lambdas[fit_succeeded] # Succeeded lambdas
  A.est <- A.est[, fit_succeeded]
  B.est <- B.est[, fit_succeeded]
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
          total_effect = sum(AB.est[, i]) + C.est[i]
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
      chosen = lambdas1 == chosen_lambda #minimum lambda meeting threshold is chosen
    )

  output$chosen_lambda = lambdas1[which_chosen]
  output$chosen_fit = results[[which_chosen]]
  output$vss = vss_out

  return(output)

}


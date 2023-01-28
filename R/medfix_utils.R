#Helper functions for implementing MedFix


get_blank <- function(X, notify = T){

  if (notify) message("   Null result. No mediators selected by penalized model.")

  out <-
    data.frame(
      mediator = colnames(X),
      beta = 0,
      sd = 0,
      p.value = 1
    )

  return(out)

}


adlasso <- function(X, Y, nlambda = 100, nlambda2 = 50, nfolds = 5){

  Y <- as.vector(Y)
  X <- as.matrix(X)
  n <- length(Y)

  if(is.null(colnames(X))) colnames(X) <- paste0("m", 1:ncol(X))

  colnames_original <- colnames(X)

  #Step 1: Obtain initial fit
  #  Choose lambda2
  if (nlambda2 > 0){
    message("   Fitting initial model with elastic net...")
    lambda2_list <- exp(seq(log(0.0001), log(0.02), length.out = nlambda2))

    nlambda2 <- length(lambda2_list)
    mincv <- numeric(nlambda2)

    for(i in 1:nlambda2){
      cv_enet <-
        cv.gcdnet(
          X,
          Y,
          lambda2 = lambda2_list[i],
          nlambda = nlambda,
          standardize = T,
          nfolds = nfolds,
          method = "ls")

      mincv[i] <- min(cv_enet$cvm)

    }

    lambda2_use <- lambda2_list[which.min(mincv)]

  }else if (nlambda2 == 0){
    message("   Fitting initial model with LASSO...")
    lambda2_use <- 0

  }else{
    stop("nlambda2 must be an integer at least 0.")

  }

  #  Fit model
  mod1 <- cv.gcdnet(X, Y, lambda2 = lambda2_use, standardize = T,
                    nfolds = nfolds, method = "ls")

  #Step 2: Compute weights and obtain second fit
  coefs <- coef(mod1, s = "lambda.min")
  if(all(coefs[-1]==0)) return(get_blank(X))
  v_star <- log(sum(coefs!=0)) / log(n)
  gamma <- ceiling(2*v_star/(1 - v_star)) + 1
  X_sd <- apply(X, 2, sd)
  weights_ad <- (abs(coefs[-1] * X_sd) + 1/n)^(-gamma)
  X1 <- X[, which(coefs[-1]!=0)]
  weights_ad_non0 <- weights_ad[which(coefs[-1] != 0)]
  message("   Fitting adaptive-weighted model...")
  mod2 <- cv.gcdnet(X1, Y, standardize = T, pf = weights_ad_non0,
                    nfolds = nfolds, method = "ls")

  #Step 3: Process results
  coefs1 <- coef(mod2, s = "lambda.min")[,1]
  if(all(coefs1[-1]==0)) return(get_blank(X))
  y_preds <- predict(mod2, X1, s = "lambda.min")
  s2 <- crossprod(Y - y_preds) / (n - sum(coefs1!=0) - 1)
  X2 <- X1[, which(coefs1[-1] != 0)]
  variance <- solve(crossprod(X2))*c(s2)
  var1 <- #effect variances
    diag(1 / apply(X2, 2, sd)) %*%
    variance %*%
    diag(1 / apply(X2, 2, sd))
  # vari <- #intercept variance
  #   apply(X2, 2, mean) %*%
  #   variance %*%
  #   apply(X2, 2, mean)

  #Empty results table
  coefs1_noint <- coefs1[-1]
  adlasso_tab <- matrix(0, sum(coefs1_noint != 0), 3)
  rownames(adlasso_tab) <- names(coefs1_noint)[which(coefs1_noint != 0)]
  colnames(adlasso_tab) <- c("beta", "sd", "p.value")

  #Fill in results
  estimates <- coefs1_noint[which(coefs1_noint != 0)]
  sds <- c(sqrt(diag(var1)))
  pvs <- 2 * (1 - pnorm(abs(estimates) / sds))
  adlasso_tab <-
    data.frame(
      variable = names(coefs1_noint[which(coefs1_noint != 0)]),
      beta = estimates,
      sd = sds,
      p.value = pvs,
      row.names = NULL
    )
  adlasso_tab$index <- which(colnames(X) %in% adlasso_tab$variable)
  #index variable used for later sorting...

  #Add in variables that were not selected
  zeros <- setdiff(colnames(X), adlasso_tab$variable)
  df_zeros <- data.frame(variable = zeros, beta = 0, sd = 0 , `p.value` = 1,
                         index = which(colnames(X) %in% zeros))
  outdat <- rbind(adlasso_tab,df_zeros)
  outdat <- outdat[order(outdat$index),] #re-order

  outdat$variable <- colnames_original
  outdat$index <- NULL

  rownames(outdat) <- NULL

  return(outdat)

}


medfix <- function(A, M, Y, C1, C2, nlambda = 100, nlambda2 = 50,
                   nfolds = 5){

  p <- ncol(M)
  if (is.null(colnames(M))) colnames(M) <- paste0("m",1:p)

  #Regress out A, C1 from Y and M
  AC1 <- data.frame(A)
  if (!is.null(C1)) AC1 <- data.frame(A, C1)
  Yr <- lm(Y ~ ., data = AC1)$residuals
  Mr <- apply(M, 2, function(x){lm(x ~ ., data = AC1)$residuals})

  #Standardize variables
  sd_Yr <- sd(Yr)
  sd_Mr <- apply(Mr, 2, sd)
  Yr1 <- as.numeric(scale(Yr))
  Mr1 <- as.matrix(scale(Mr))

  #Fit adaptive LASSO on the outcome model
  message("Fitting outcome model with adaptive LASSO...")
  ym_out <- adlasso(Mr1, Yr1, nlambda = nlambda, nlambda2 = nlambda2,
                     nfolds = nfolds)[,c(1,2,4)]
  colnames(ym_out) <- c("mediator","beta","beta_pv")

  #Fit mediator models
  message("Fitting mediator models...")
  AC2 <- data.frame(A)
  if (!is.null(C2)) AC2 <- data.frame(A, C2)
  ma_out <-
    lapply(as.data.frame(M),
           function(x){summary(lm(x ~ ., data = AC2))$coefficients[2,c(1,4)]}
    )
  ma_out <- as.data.frame(do.call(rbind, ma_out))
  colnames(ma_out) <- c("alpha","alpha_pv")

  #Summarize output
  out <-
    data.frame(
      mediator = ym_out$mediator,
      beta = ym_out$beta * sd_Yr / sd_Mr,
      alpha = ma_out$alpha,
      alpha_beta = ma_out$alpha * ym_out$beta * sd_Yr / sd_Mr,
      ab_pv = pmax(ym_out$beta_pv,ma_out$alpha_pv),
      te = summary(lm(Y ~ ., data = AC1))$coefficients[2,1],
      row.names = NULL
    )

  return(out)


}


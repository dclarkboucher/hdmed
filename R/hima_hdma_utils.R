#Helper functions for hdma.R and hima.R

#Sources:
#https://github.com/YinanZheng/HIMA
#https://github.com/YuzhaoGao/High-dimensional-mediation-analysis-R

blank_effects <- function(){
  warning("Null result. No mediators selected by the penalized model...")
  data.frame(matrix(NA,0,6))
}

hima <- function(X, Y, M, COV.XM = NULL, COV.MY = COV.XM,
                 Y.family = c("gaussian", "binomial"),
                 M.family = c("gaussian", "negbin"),
                 penalty = c("MCP", "SCAD", "lasso"),
                 topN = NULL,
                 verbose = FALSE,
                 ...) {

  Y.family <- match.arg(Y.family)
  M.family <- match.arg(M.family)
  penalty <- match.arg(penalty)

  n <- nrow(M)
  p <- ncol(M)

  if(is.null(topN)) {
    if(Y.family == "binomial") d <- ceiling(n/(2*log(n))) else d <- ceiling(2 * n/log(n))
  } else {
    d <- topN  # the number of top mediators that associated with exposure (X)
  }

  d <- min(p, d) # if d > p select all mediators

  #Step 1: SIS
  message("Screening mediators...")
  if(Y.family == "binomial")
  {
    # Screen M using X given the limited information provided by Y (binary)
    if(verbose) message("    Screening M using the association between X (independent variable) and M (dependent variable): ", appendLF = FALSE)
    alpha = SIS_Results <- himasis(NA, M, X, COV.XM, glm.family = M.family, modelstatement = "Mone ~ X",
                                   verbose, tag = paste0("Sure Independent Screening (M ~ X + COV.XM, family: ", M.family, ")"))
    SIS_Pvalue <- SIS_Results[2,]
  } else if(Y.family == "gaussian"){
    # Screen M using Y (continuous)
    if(verbose) message("    Screening M using the association between M (independent variable) and Y (dependent variable): ", appendLF = FALSE)
    SIS_Results <- himasis(Y, M, X, COV.MY, glm.family = Y.family, modelstatement = "Y ~ Mone + X",
                           verbose, tag = paste0("Sure Independent Screening (Y ~ M + X + COV.MY, family: ", Y.family, ")"))
    SIS_Pvalue <- SIS_Results[2,]
  } else {
    stop(paste0("Family ", Y.family, " is not supported."))
  }

  # Note: ranking using p on un-standardized data is equivalent to ranking using beta on standardized data
  SIS_Pvalue_sort <- sort(SIS_Pvalue)
  ID <- which(SIS_Pvalue <= SIS_Pvalue_sort[d])  # the index of top mediators
  if(verbose) message("    Top ", length(ID), " mediators are selected: ", paste0(names(SIS_Pvalue_sort[seq_len(d)]), collapse = ","))

  M_SIS <- M[, ID]
  XM <- cbind(M_SIS, X)

  #Step 2: MCP
  message("Fitting outcome model with MCP...")
  ## Based on the screening results in step 1. We will find the most influential M on Y.
  if(is.null(COV.MY)) {
    fit <- ncvreg(XM, Y, family = Y.family,
                  penalty = penalty,
                  penalty.factor = c(rep(1, ncol(M_SIS)), 0), ...)
  } else {
    COV.MY <- data.frame(COV.MY)
    COV.MY <- data.frame(model.matrix(~., COV.MY))[, -1]
    conf.names <- colnames(COV.MY)
    if(verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    XM_COV <- cbind(XM, COV.MY)
    fit <- ncvreg(XM_COV, Y, family = Y.family,
                  penalty = penalty,
                  penalty.factor = c(rep(1, ncol(M_SIS)), rep(0, 1 + ncol(COV.MY))), ...)
  }
  # plot(fit)

  lam <- fit$lambda[which.min(BIC(fit))]
  if(verbose) message("    Tuning parameter lambda selected: ", lam)
  Coefficients <- stats::coef(fit, lambda = lam)
  est <- Coefficients[2:(d + 1)]
  ID_1_non <- which(est != 0)
  if(length(ID_1_non) == 0) return(blank_effects())

  #Step 3: Mediator models
  message("Fitting mediator models...")
  if(verbose) message("    Non-zero ", penalty, " beta estimate(s) of mediator(s) found: ", paste0(names(ID_1_non), collapse = ","))
  beta_est <- est[ID_1_non]  # The non-zero MCP estimators of beta
  ID_test <- ID[ID_1_non]  # The index of the ID of non-zero beta in Y ~ M
  ##

  if(Y.family == "binomial")
  {
    ## This has been done in step 1 (when Y is binary, alpha is estimated in M ~ X)
    alpha <- alpha[, ID_test, drop = FALSE]
    message("    Using alpha estimated in Step 1 ...   (", Sys.time(), ")")
  } else if(Y.family == "gaussian"){
    if(verbose) message("    Estimating alpha (effect of X on M): ", appendLF = FALSE)
    alpha <- himasis(NA, M[, ID_test, drop = FALSE], X, COV.XM, glm.family = M.family,
                     modelstatement = "Mone ~ X", verbose,
                     tag = paste0("site-by-site alpha estimation (M ~ X + COV.XM, family: ", M.family, ")"))
  } else {
    stop(paste0("Family ", Y.family, " is not supported."))
  }


  alpha_est_ID_test <- as.numeric(alpha[1, ])  #  the estimator for alpha
  P_alpha <- alpha[2, ]  # the adjusted p-value for alpha (bonferroni)

  alpha_est <- alpha_est_ID_test

  ## Post-test based on the oracle property of the MCP penalty
  if(is.null(COV.MY)) {
    YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X)
  } else {
    YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X, COV.MY)
  }

  res <- summary(glm(Y ~ ., family = Y.family, data = YMX))$coefficients
  est <- res[2:(length(ID_test) + 1), 1]  # the estimator for beta
  P_beta <- res[2:(length(ID_test) + 1), 4]
  ab_est <- alpha_est * beta_est

  ## Use the maximum value as p value
  P_value <- pmax(P_beta, P_alpha)

  # Total effect
  if(is.null(COV.MY)) {
    YX <- data.frame(Y = Y, X = X)
  } else {
    YX <- data.frame(Y = Y, X = X, COV.MY)
  }

  gamma_est <- stats::coef(glm(Y ~ ., family = Y.family, data = YX))[2]

  results <-
    data.frame(
      mediator = colnames(M_SIS)[ID_1_non],
      alpha = alpha_est,
      alpha_pv = P_alpha,
      beta = beta_est,
      beta_pv = P_beta,
      alpha_beta = ab_est,
      ab_pv = P_value,
      te = gamma_est,
      row.names = NULL
    )

  return(results)

}

hdma <- function (X, Y, M, COV.XM = NULL, COV.MY = COV.XM,
                  family = c("gaussian","binomial"),
                  method = c("lasso", "ridge"), topN = NULL,
                  verbose = F, ...){

  family <- match.arg(family)
  method <- match.arg(method)
  n <- nrow(M)
  p <- ncol(M)
  if (is.null(topN)) {
    if (family == "binomial") d <- ceiling(n/(2*log(n))) else d <- ceiling(2*n/log(n))
  } else {
    d <- topN      # the number of top mediators that associated with independent variable (X)
  }
  d <- min(p, d)   # if d > p select all mediators


  #Step 1: SIS
  message("Screening mediators...")
  if(family == "binomial")
  {
    alpha = SIS_Results <- himasis(NA, M, X, COV.XM, glm.family = "gaussian", modelstatement = "Mone ~ X", verbose, tag = "Sure Independent Screening")
    SIS_Pvalue <- SIS_Results[2,]
  } else if (family == "gaussian"){
    # Screen M using Y (continuous)
    if(verbose) message("Screening M using the association between M and Y: ", appendLF = FALSE)
    SIS_Results <- himasis(Y, M, X, COV.MY, glm.family = family, modelstatement = "Y ~ Mone + X", verbose, tag = "Sure Independent Screening")
    SIS_Pvalue <- SIS_Results[2,]
  } else {
    stop(paste0("Family ", family, " is not supported."))
  }
  # Note: ranking using p on un-standardized data is equivalent to ranking using beta on standardized data
  SIS_Pvalue_sort <- sort(SIS_Pvalue)
  ID <- which(SIS_Pvalue <= SIS_Pvalue_sort[d])  # the index of top mediators
  if(verbose) message("Top ", length(ID), " mediators are selected: ", paste0(names(SIS_Pvalue_sort[seq_len(d)]), collapse = ","))
  M_SIS <- M[, ID]
  XM <- cbind(M_SIS, X)

  #Step 2: LASSO
  message("Fitting outcome model with de-biased LASSO...")
  ## Based on the SIS results in step 1. We will find the most influential M on Y.
  if (is.null(COV.MY)) {
    set.seed(1029)
    if (method == "lasso") fit <- lasso.proj(XM, Y, family = family) else fit <- ridge.proj(XM, Y, family = family)
  } else {
    COV.MY <- data.frame(COV.MY)
    COV.MY <- data.frame(model.matrix(~., COV.MY))[, -1]
    conf.names <- colnames(COV.MY)
    XM_COV <- cbind(XM, COV.MY)
    set.seed(1029)
    if (method == "lasso") fit <- lasso.proj(XM_COV, Y, family = family) else fit <- ridge.proj(XM_COV, Y, family = family)
  }
  P_hdi <- fit$pval[1:length(ID)]
  index <- which(P_hdi<=0.05)

  if(length(index)==0) return(blank_effects())

  ID_test <- ID[index]

  #Step 3: Mediator models
  message("Fitting mediator models...")
  if(family == "binomial")
  {
    ## This has been done in step 1 (when Y is binary)
    alpha_est <- alpha[,ID_test, drop = FALSE]
  } else {
    alpha_est <- himasis(NA, M[,ID_test, drop = FALSE], X, COV.XM, glm.family = "gaussian", modelstatement = "Mone ~ X",
                         verbose, tag = "site-by-site ordinary least squares estimation")
  }
  beta_P<-P_hdi[index]
  beta_hat<-fit$bhat[index]              # the estimator for beta
  alpha_hat<-as.numeric(alpha_est[1, ])
  ab_est<-beta_hat*alpha_hat
  alpha_P<-alpha_est[2,]
  PA <- rbind(beta_P,alpha_P)
  P_value <- apply(PA,2,max)

  if (is.null(COV.MY)) {
    YX <- data.frame(Y = Y, X = X)
  } else {
    YX <- data.frame(Y = Y, X = X, COV.MY)
  }
  gamma_est <- stats::coef(glm(Y ~ ., family = family, data = YX))[2]
  results <-
    data.frame(
      mediator = colnames(M_SIS)[index],
      alpha = alpha_hat,
      alpha_pv = alpha_P,
      beta = beta_hat,
      beta_pv = beta_P,
      alpha_beta = ab_est,
      ab_pv = P_value,
      te = gamma_est,
      row.names = NULL
    )

  return(results)
}

#Function for sure independent screening
globalVariables("n")
globalVariables("M_chunk")

himasis <- function(Y, M, X, COV, glm.family, modelstatement, verbose, tag) {
  L.M <- ncol(M)
  M.names <- colnames(M)

  X <- data.frame(X)
  X <- data.frame(model.matrix(~., X))[, -1]

  if (is.null(COV)) {
    if (verbose) message("    No covariate is adjusted")
    datarun <- data.frame(Y = Y, Mone = NA, X = X)
    modelstatement <- modelstatement
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~., COV))[, -1]
    conf.names <- colnames(COV)
    if (verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    datarun <- data.frame(Y = Y, Mone = NA, X = X, COV = COV)
    modelstatement <- eval(parse(text = (paste0(modelstatement, "+",
                                                paste0(paste0("COV.", conf.names), collapse = "+")))))
  }

  doOne <- eval(parse(text = doOneGen(paste0("try(glm(modelstatement, family = ",
                                             glm.family, ", data = datarun))"), "c(1,4)")))

  registerDoSEQ()
  ncore <- 1
  results <- foreach(n = iterators::idiv(L.M, chunks = ncore),
                     M_chunk = iblkcol_lag(M, chunks = ncore),
                     .combine = "cbind") %dopar% {sapply(seq_len(n), doOne, datarun, M_chunk)}

  colnames(results) <- M.names
  return(results)
}

#Helper function used by himasis (above)
doOneGen <- function(model.text, colind.text) {
  L <- length(eval(parse(text = colind.text)))
  script <- paste0("doOne <- function(i, datarun, Ydat){datarun$Mone <- Ydat[,i]; model <- ",
                   model.text, ";if('try-error' %in% class(model)) b <- rep(NA, ",
                   L, ") else { res=summary(model)$coefficients; b <- res[2,", colind.text,
                   "]};invisible(b)}")
  return(script)
}

#Helper function used by himasis (above)
iblkcol_lag <- function(M, ...) {
  i <- 1
  it <- iterators::idiv(ncol(M), ...)

  nextEl <- function() {
    n <- iterators::nextElem(it)
    r <- seq(i, length = n)
    i <<- i + n
    M[, r, drop = FALSE]
  }
  obj <- list(nextElem = nextEl)
  class(obj) <- c("abstractiter", "iter")
  obj
}



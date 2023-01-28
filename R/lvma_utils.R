# Utility functions for LVMA


# LVMA function called by mediate_lvma wrapper
lvma <- function (A, M, Y, q, rhoL, rhoE, rhoY, imax) {

      Obs_Data = list(
        E = matrix(A, ncol = 1),
        Y = matrix(Y, ncol = 1),
        M = M,
        Covariates = NULL
      )

      MR <- apply(M, 2, function(x) lm(x ~ A)$residuals)
      aa = prcomp(MR, scale = TRUE, center = TRUE)
      Loadings = varimax(aa$rotation[, 1:q])$loadings

      out <-
        run_LVMA_Cont(
          Obs_Data = Obs_Data,
          Loadings,
          rhoL = rhoL,
          rhoE = rhoE,
          rhoY = rhoY,
          niter = imax
        )
  out

}


process_lvma <- function (out, type = c("EBIC", "BIC", "AIC"), M, rhoLM, rhoEL, rhoLY){

  type <- match.arg(type)
  value <- out[[paste0("val_", type)]]

  # Which tuning parameters were chosen?
  which_chosen <- out[[paste0("model_selection",type)]]
  penalty <- c(rhoLM[which_chosen[1]], rhoEL[which_chosen[2]],
               rhoLY[which_chosen[3]])
  names(penalty) <- c("rhoLM", "rhoEL", "rhoLY")

  # What are the estimated effects?
  solution <- out[[paste0(type,"_sol")]]
  LM_effects <- solution$Lambda

  q <- ncol(LM_effects)
  colnames(LM_effects) <- paste0("L", 1:q)
  rownames(LM_effects) <- colnames(M)

  LY_effects <- as.vector(solution$betaFY)
  names(LY_effects) <- paste0("L", 1:q)

  AL_effects <- as.vector(solution$betaEF)
  names(AL_effects) <- paste0("L", 1:q)

  AY_direct_effect <- solution$betaEY

  # Which mediators are active?
  mediator_active <-
    suppressWarnings(apply(t(LM_effects) * (AL_effects * LY_effects), 2, any)) |>
    as.numeric()
  names(mediator_active) <- colnames(M)

  output <-
    list(
      penalty = penalty,
      AL_effects = AL_effects,
      LM_effects = as.data.frame(LM_effects),
      AY_direct_effect = AY_direct_effect,
      LY_effects = LY_effects,
      mediator_active
    )

  output[[type]] <- value

  return(output)
}



# Functions for LVMA - continuous -----------------------------------------

run_LVMA_Cont <- function(Obs_Data = Obs_Data,
                          Loadings,
                          rhoL,
                          rhoE,
                          rhoY,
                          niter = 5000) {
  E = matrix(Obs_Data$E, ncol = 1)
  M = Obs_Data$M
  Cova = Obs_Data$Covariates
  Y = matrix(Obs_Data$Y, ncol = 1)
  M = Obs_Data$M
  nf = ncol(Loadings)

  betaEY = 0
  sigmay = var(Y)
  betaFE = matrix(rep(0, nf), ncol = 1)
  betaFY = matrix(rep(0, nf), ncol = 1)
  MR = NULL
  for (i in 1:length(M[1, ])) {
    MR = cbind(MR, lm(M[, i] ~ E)$residuals)
  }
  nf = length(Loadings[1, ])
  Lambda = Loadings
  psiM = diag(apply(M, 2, var))
  Means = colMeans(M)
  meansY = mean(Y)
  rho10 = Lambda * 0 + 0
  rho20 = (betaFE) * 0 + 0
  rho30 = c(t(betaFE) * 0 + 0)

  if (is.null(Cova)) {
    betaCovF = NULL
    betaCovY = NULL
  } else{
    betaCovF = (matrix(0.01, ncol = length(Cova[1, ]), nrow = length(Lambda[1, ])))
    betaCovY = rep(0.01, length(Cova[1, ]))
  }


  Fit0 = get_estALS_Cont(
    E,
    M,
    Y,
    betaFE,
    betaEY,
    betaFY,
    Lambda,
    psiM,
    sigmay,
    Means,
    meansY,
    rho10,
    rho20,
    rho30,
    niter,
    Cova,
    betaCovF,
    betaCovY
  )

  betaEF = betaEF0 = Fit0$betaEF
  betaFY = betaFY0 = Fit0$betaFY
  betaEY = Fit0$betaEY
  Lambda = Lambda0 = Fit0$Lambda
  sigmaY = Fit0$sigmaY
  psiM = Fit0$ThetaM
  meansY = Fit0$muY
  Means = Fit0$meansM
  betaCovF = Fit0$betaCovF
  betaCovX = Fit0$betaCovX

  BB = fit_ALS_Cont(
    E,
    M,
    Y,
    betaEF,
    betaEY,
    betaFY,
    Lambda,
    psiM,
    sigmay,
    Means,
    meansY,
    rhoL,
    rhoE,
    rhoY,
    Lambda0,
    betaEF0,
    betaFY0,
    niter,
    Cova,
    betaCovF,
    betaCovY
  )
  return(BB)
}

fit_ALS_Cont <- function(E,
                         M,
                         Y,
                         betaEF,
                         betaEY,
                         betaFY,
                         Lambda,
                         psiM,
                         sigmay,
                         Means,
                         meansY,
                         rhoL,
                         rhoE,
                         rhoY,
                         Lambda0,
                         betaEF0,
                         betaFY0,
                         niter = 5000,
                         Covariates = NULL,
                         betaCovF = NULL,
                         betaCovY = NULL) {
  BIC_sol = NULL
  EBIC_sol = NULL
  AIC_sol = NULL

  valBIC = 10 ^ 10
  valEBIC = 10 ^ 10
  valAIC =  10 ^ 10

  selection_BIC = NULL
  selection_EBIC = NULL
  selection_AIC = NULL

  iL = length(rhoL)
  iE = length(rhoE)
  iY = length(rhoY)

  StartI = StartBetaY = StartLambda = StartBetaE = list(
    meansM = Means,
    muY = meansY,
    betaEY = betaEY,
    betaFY = betaFY,
    Lambda = Lambda,
    betaEF = betaEF,
    sigmaY = sigmay,
    ThetaM = psiM,
    betaCovF = betaCovF,
    betaCovY = betaCovY
  )


  for (i in 1:iL) {
    StartI = StartLambda
    for (j in 1:iE) {
      StartI = StartBetaE
      for (k in 1:iY) {

        AA = get_estALS_Cont(
          E,
          M,
          Y,
          StartI$betaEF,
          StartI$betaEY,
          StartI$betaFY,
          StartI$Lambda,
          StartI$ThetaM,
          StartI$sigmaY,
          StartI$meansM,
          StartI$muY,
          rhoL[i] / abs(Lambda0),
          (rhoE[j] / abs(betaEF0)),
          matrix(c(rhoY[k] / abs(betaFY0)), nrow = 1),
          niter,
          Covariates,
          StartI$betaCovF,
          StartI$betaCovY
        )
        if (k == 1) {
          StartBetaE = AA
        }

        StartI = AA

        X = cbind(Covariates, E)
        Val = get_Penalty(X, M, Y, AA)

        if (Val$BIC <= valBIC) {
          BIC_sol = AA
          valBIC = Val$BIC
          selection_BIC = c(rhoL[i], rhoE[j], rhoY[k])
        }
        if (Val$EBIC <= valEBIC) {
          EBIC_sol = AA
          valEBIC = Val$EBIC
          selection_EBIC = c(rhoL[i], rhoE[j], rhoY[k])
        }
        if (Val$AIC <= valAIC) {
          AIC_sol = AA
          valAIC = Val$AIC
          selection_AIC = c(rhoL[i], rhoE[j], rhoY[k])
        }
        if (sum(AA$betaFY != 0) == 0) {
          break
        }
      }
      if (sum(AA$betaEF != 0) == 0) {
        break
      }
    }
  }
  return (
    list(
      model_selectionEBIC = selection_EBIC,
      EBIC_sol = EBIC_sol,
      val_EBIC = valEBIC,
      model_selectionBIC = selection_BIC,
      BIC_sol = BIC_sol,
      val_BIC = valBIC,
      model_selectionAIC = selection_AIC,
      AIC_sol = AIC_sol,
      val_AIC = valAIC
    )
  )
}

get_estALS_Cont <- function(E,
                            M,
                            Y,
                            betaEF,
                            betaEY,
                            betaFY,
                            Lambda,
                            psiM,
                            sigmay,
                            Means,
                            meansY,
                            rho1,
                            rho2,
                            rho3,
                            niter = 5000,
                            Covariates = NULL,
                            betaCovF = NULL,
                            betaCovY = NULL) {
  rho2 = matrix(rho2, ncol = 1)
  rho3 = matrix(c(0, rho3), nrow = 1)

  ##
  ## Combine Exposure with covariates
  ##
  E = cbind(Covariates, E)
  if (!is.null(Covariates)) {
    rho2 = cbind(matrix(
      0,
      nrow = length(betaCovF[, 1]),
      ncol = length(Covariates[1, ])
    ), rho2)
    rho3 = cbind(rep(0, length(betaCovY)), rho3)
  }

  # nun of cov. + 1 by num of factors   matrix of effects on factors
  betaEF = cbind(betaCovF, (betaEF))
  betaY = c(betaCovY, betaEY, betaFY)
  betaEY = betaY[(1:length(E[1, ]))]

  #
  # Initiation
  #


  theta = c(betaY,
            as.vector(betaEF),
            sigmay,
            as.vector(psiM),
            as.vector(Lambda)) # initial value
  thetaN = theta + 0.1
  betaFEEst = betaEF
  betaYEst = betaY
  betaEYEst = betaEY
  betaFYEst = betaFY
  lambdaEst = Lambda
  sigmayEst = sigmay
  psiMEst = psiM

  mu0MEst = Means
  mu0YEst = meansY

  kkkIt = 1
  Val1 = 1
  ccc = 1

  # message('Maximizing joint likelihood')


  while (max(abs(theta - thetaN)) > 0.0000001 & (ccc) > 0.001) {
    MeansMat = get_MEANS(M, Y, E, lambdaEst, betaFEEst, betaYEst, mu0MEst, mu0YEst)
    MeansY = MeansMat$muY
    MeansF = MeansMat$muF
    Means = MeansMat$muM

    #
    # Set-up covarinces
    #
    SigmaFF = diag(1, length(lambdaEst[1, ]))
    VarM = lambdaEst %*% t(lambdaEst)  + psiMEst
    CovMY = lambdaEst %*% betaFYEst
    CovYM = t(CovMY)
    varY = t(betaFYEst) %*% betaFYEst + sigmayEst
    covFM = t(lambdaEst)
    covFY = betaFYEst
    SigmaEstMY = rbind(cbind(VarM, CovMY), cbind(CovYM, varY))
    SigmaEstFMY = cbind(covFM, covFY)

    FF = getF(M, Y, MeansF, Means, MeansY, SigmaEstFMY, SigmaEstMY, SigmaFF)### Calculate E(Factor|DATA) and Calculate E(Factor'Factor|DATA)

    Est0 = regMonFALS(M, FF, rho1, psiMEst, lambdaEst, mu0MEst) # estimate effects of F on M
    lambdaEst = Est0$lambda
    psiMEst = Est0$sigmaM
    mu0MEst = Est0$means

    Est1 = regFonEALS(E, FF, (rho2), betaFEEst) # estimate effects of F on E
    betaFEEst = (Est1$rho)
    Est2 = regFEonYALS(Y, E, FF, rho3, sigmayEst, betaYEst, mu0YEst) # estimate effects of F and E on Y
    mu0YEst = Est2$mu0
    betaYEst = Est2$beta
    sigmayEst = Est2$sigma


    thetaN = theta
    theta = c(
      betaYEst,
      as.vector(betaFEEst),
      sigmayEst,
      as.vector(psiMEst),
      as.vector(lambdaEst)
    )

    kkkIt = kkkIt + 1
    if (kkkIt >= niter) {
      break
    }

    if (is.null(betaCovY)) {
      AA = list(
        meansM = mu0MEst,
        muY = mu0YEst,
        betaEY = betaYEst[1],
        betaFY = betaYEst[-1],
        Lambda = lambdaEst,
        betaEF = betaFEEst,
        sigmaY = sigmayEst,
        ThetaM = psiMEst,
        betaCovF = NULL,
        betaCovY = NULL
      )
    } else{
      kkk = 1:length(E[1, ])
      betaFMatr = matrix(betaFEEst[, -length(kkk)], ncol = length(kkk) -
                           1)
      betaFEEstM = matrix(betaFEEst[, length(kkk)], ncol = 1)
      AA = list(
        meansM = mu0MEst,
        muY = mu0YEst,
        betaEY = betaYEst[kkk][length(kkk)],
        betaFY = betaYEst[-kkk],
        Lambda = lambdaEst,
        betaEF = betaFEEstM,
        sigmaY = sigmayEst,
        ThetaM = psiMEst,
        betaCovF = betaFMatr,
        betaCovY = betaYEst[kkk][-length(kkk)]
      )
    }

    Val = get_Penalty(E, M, Y, AA)

    Val0 = Val1

    Val1 = Val$logL + sum(rho1 * abs(lambdaEst)) + sum(as.vector(rho2) *
                                                         abs(betaFEEst)) + sum(as.vector(rho3) * abs(betaYEst[]))
    ccc = 1
    if (kkkIt > 10) {
      ccc = Val0 - Val1
    }

  }
  # message('Done maximizing')
  return (AA)
}


get_Penalty <- function(E, M, Y, AA) {
  N = length(Y)
  lambdaEst = Lambda = AA$Lambda
  betaEYst = betaEY = AA$betaEY
  betaEFEst = betaFE = (cbind(AA$betaCovF, AA$betaEF))

  betaFYEst = betaFY = AA$betaFY
  betaYst  = betaY = c(AA$betaCovY, betaEYst, AA$betaFY)
  mu0MEst = AA$meansM
  mu0YEst = AA$muY
  sigmayEst = AA$sigmaY
  psiMEst = AA$ThetaM

  MeansMat = get_MEANS(M, Y, E, lambdaEst, betaEFEst, betaYst, mu0MEst, mu0YEst)
  MeansY = MeansMat$muY
  MeansF = MeansMat$muF
  Means = MeansMat$muM

  MU  = cbind(Means, MeansY)
  VarM = lambdaEst %*% t(lambdaEst)  + psiMEst
  CovMY = lambdaEst %*% betaFYEst
  CovYM = t(CovMY)
  varY = betaFYEst %*% betaFYEst + sigmayEst
  covFM = t(lambdaEst)
  covFY = betaFYEst

  SigmaEstMY = rbind(cbind(VarM, CovMY), cbind(CovYM, varY))
  Data = cbind(M, Y)

  DATA = Data  - MU
  DataCase = as.matrix(DATA)
  Scase = t(DataCase) %*% DataCase / N

  LLcase = N / 2 * (log(det(SigmaEstMY)) + tr(solve(SigmaEstMY) %*% Scase))
  LL = LLcase
  llAIC = 2 * LL + 2 * (sum(betaFE != 0) + sum(betaFY != 0) + sum(betaEY !=
                                                                    0) + sum(Lambda != 0))
  llBIC = 2 * LL + (sum(betaFE != 0) + sum(betaFY != 0) + sum(betaEY != 0) + sum(Lambda !=
                                                                                   0)) * log(N)
  bb1 = sum(Lambda != 0) + sum(betaFY != 0) + sum(betaFE != 0)
  nn1 = length(Lambda != 0) + length(betaFY != 0) + length(betaFE != 0)
  BB1 = sum(log((nn1 - bb1 + 1):nn1)) - sum(log(1:bb1))

  llEBIC = 2 * LL + (sum(betaFE != 0) + sum(betaFY != 0) + sum(betaEY !=
                                                                 0) + sum(Lambda != 0)) * log(N) + BB1
  return(list(
    logL = LL,
    AIC = llAIC,
    BIC = llBIC,
    EBIC = llEBIC
  ))

}

get_MEANS <- function(M,
                      Y,
                      E,
                      Lambda,
                      betaEF,
                      betaEY,
                      Means,
                      MeansY) {
  L = length(Y)
  muF = E %*% t(betaEF)

  CommonM = t(matrix(Means, nrow = length(Means), ncol = L))
  muM  = CommonM + muF %*% t(Lambda)

  CommonY = rep(MeansY, L)
  muY = CommonY + cbind(E, muF) %*% betaEY

  return (list(muF = muF, muY = muY, muM = muM))
}

getF <-
  function(M,
           Y,
           MeansF,
           Means,
           MeansY,
           SigmaEstFMY,
           SigmaEstMY,
           SigmaFF) {
    Data = rbind(t(M - Means), t(Y - MeansY))
    F = t(MeansF) + SigmaEstFMY %*% solve(SigmaEstMY) %*% Data
    VarFF = SigmaFF - SigmaEstFMY %*% solve(SigmaEstMY) %*% t(SigmaEstFMY)
    FF = VarFF * length(Y) + F %*% t(F)
    return (list(f = F, ff = FF))
  }

regMonFALS = function(M, FF, rho, psiMEst, Lambda, Means) {
  #solve for lambda
  ThetaM = diag(psiMEst)
  Y = M
  X = FF$f
  L = length(X[1, ])
  X = rbind(rep(1, L), X)
  #lambda = matrix(t(Inv%*%X%*%Y),ncol=length(X[,1]))
  n = length(FF$f[1, ])
  xca = as.vector(rowSums(FF$f))
  AA = rbind(xca, FF$ff)
  AA = cbind(c(n, xca), AA)
  rho = cbind(rep(0, length(rho[, 1])), rho)

  lambda = t(get_regressionALS(AA, X, Y, ThetaM, rho, t(cbind(Means, Lambda))))
  ThetaM = NULL
  means = lambda[, 1]
  lambda = lambda[, -1]
  X = X[-1, ]
  for (i in 1:length(lambda[, 1])) {
    Y[, i] = Y[, i] - means[i]

    ThetaM = c(ThetaM, (
      sum(Y[, i] ^ 2) + t(lambda[i, ]) %*% (FF$ff) %*% (lambda[i, ]) - 2 * t(Y[, i]) %*%
        t(X) %*% (lambda[i, ])
    ) / length(Y[, 1]))
  }
  #save(rho,FF,lambda,means,ThetaM,psiMEst,file='rho')
  return (list(
    lambda = lambda,
    sigmaM = diag(ThetaM),
    means = means
  ))
}


regFonEALS = function(X, FF, rho, Beta) {
  X = matrix(X, ncol = length(X[1,]))
  F = FF$f
  Y = t(cbind(F))
  sigmaF = rep(1, length(Y[1,]))
  rhoe = get_regressionALS(t(X) %*% X, t(X), Y, sigmaF, (rho), t(Beta))
  return (list(rho = t(rhoe)))
}

regFEonYALS = function(Y, E, FF, rho, sigmaY, Beta, meansY) {
  Y = matrix(Y, ncol = 1)
  E = matrix(E, ncol = length(E[1, ]))
  X = FF$f
  L = length(X[1, ])
  X = rbind(rep(1, L), t(E), X)
  #lambda = matrix(t(Inv%*%X%*%Y),ncol=length(X[,1]))
  n = length(FF$f[1, ])
  xca = as.vector(rowSums(FF$f))

  AA = FF$ff
  BB = FF$f %*% E
  AA = rbind(xca, t(BB), AA)
  CC = rbind(c(n, t(t(E) %*% rep(1, L))), cbind(t(E) %*% rep(1, L), t(E) %*%
                                                  E))
  AA = cbind(rbind(CC, cbind(xca, BB)), AA)
  rho = matrix(c(0, rho), nrow = 1)

  Beta = matrix(c(meansY, as.vector(Beta)), ncol = 1)
  betaY = (get_regressionALS(AA, X, Y, sigmaY, rho, Beta))

  sigma = (sum(Y ^ 2) + t(betaY) %*% AA %*% (betaY) - 2 * t(Y) %*% t(X) %*%
             betaY) / length(Y[, 1])

  return (list(
    mu0 = betaY[1],
    beta = betaY[-1],
    sigma = sigma
  ))
}

get_regressionALS = function(XtX, Xt, Y, Sigma, rho, BetaM) {
  numY = length(Y[1, ]) #
  Beta = NULL
  for (i in 1:numY) {
    beta = get_single_regressionALS(XtX, Xt, Y[, i], Sigma[i], rho[i, ], BetaM[, i])
    Beta = cbind(Beta, beta)
  }
  Beta[abs(Beta) < 10 ^ (-16)] = 0
  return (Beta)
}

get_single_regressionALS = function(XtX, Xt, Y, Sigma, rho, Beta) {
  L = length(XtX[1, ])

  beta0 = matrix(c(rep(0, L - length(Beta)), as.vector(Beta)), ncol = 1)
  beta1 = beta0 + 0.01
  kkk = 0
  XtY =  Xt %*% Y
  while (sum(abs(beta0 - beta1)) > 0.00000001) {
    beta1 = beta0
    kkk = kkk + 1
    for (i in 1:L) {
      est = solve(XtX[i, i]) %*% (XtY[i] - XtX[i, -i] %*% beta0[-i, ])
      rhostr = rho[i] * Sigma / XtX[i, i]
      beta0[i, ] = soft_thrALS(rhostr, est)
    }
  }
  return (beta0)
}

soft_thrALS = function(rhostr, est) {
  if ((abs(est) - rhostr) < 0) {
    return (0)
  }

  x = sign(est) * (abs(est) - rhostr)
  return (x)
}

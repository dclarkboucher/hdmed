
pca_loadings <-
  function(X, M, adaptive = FALSE, var.per = 0.8, n.pc=NULL)
  {

    n <- nrow(M)
    p <- ncol(M)
    if (is.null(colnames(M))) colnames(M) <- paste0("m",1:p)

    # PCA on M~X residuals
    fit.m <- lm(M ~ X)
    E.m <- fit.m$residuals
    Sigma.m <- cov(E.m)
    svd.m <- svd(Sigma.m)

    #How many PCs?
    if (adaptive){
      n.pc <- sum(cumsum(svd.m$d) / sum(svd.m$d) < var.per) + 1

    }else{
      if (is.null(n.pc)){
        n.pc <- min(p,n)
      }
    }

    #Loading matrix U
    U <- svd.m$u[, 1:n.pc]
    colnames(U) <- paste0("u", 1:n.pc)
    rownames(U) <- colnames(M)

    var.pc <- cumsum(svd.m$d[1:n.pc]) / sum(svd.m$d)
    out <- list(U = U, var.pc = var.pc)
    return(out)
  }

spca_loadings <- #Fixed!
  function(X, M, adaptive = FALSE, fused = FALSE, var.per = 0.8, n.pc=NULL,
           gamma = 0, eps = 1e-4, trace = TRUE, maxsteps = 2000,
           lambda.tune = c("R2"), per.jump = 0.7)
  {

    n <- nrow(M)
    p <- ncol(M)
    if (is.null(colnames(M))) colnames(M) <- paste0("m",1:p)

    # PCA on M~X residuals
    fit.m <- lm(M ~ X)
    E.m <- fit.m$residuals
    Sigma.m <- cov(E.m)
    svd.m <- svd(Sigma.m)

    #How many PCs?
    if (adaptive){
      n.pc <- sum(cumsum(svd.m$d) / sum(svd.m$d) < var.per) + 1

    }else{
      if (is.null(n.pc)){
        n.pc <- min(p,n)
      }
    }

    U <- svd.m$u[, 1:n.pc]
    colnames(U) <- paste0("PC", 1:n.pc)
    rownames(U) <- colnames(M)

    var.pc <- cumsum(svd.m$d[1:n.pc]) / sum(svd.m$d)

    # Sparse PCA
    if (lambda.tune[1]=="R2"){

      #Check for fused LASSO
      if (fused){
        D <- NULL
      }else{
        D <- diag(1,ncol(E.m))
      }

      #Chose lambda by adjusted total variance
      re.SPCA <- R2.flasso(E.m, U, D = D, gamma = gamma, eps=eps,
                           maxsteps = maxsteps, per.jump = per.jump)
      V <- re.SPCA$V
      SPCA.var.per.cum <- re.SPCA$var.per
    }

    W <- apply(V, 2, function(x){return(x / sqrt(sum(x^2)))})
    colnames(V) = colnames(W) <- paste0("PC",1:n.pc)
    rownames(V) = rownames(W) <- colnames(M)

    # U: original loading matrix
    # V: sparsify loading matrix
    # W: sparsify loading matrix with l2-norm 1
    out <- list(U = U, V = V, W = W, var.pc = var.pc,
                var.spc = SPCA.var.per.cum)
    return(out)
  }


R2.flasso <- #Fixed!
  function(E,U,D=NULL,gamma=0,eps=1e-4,maxsteps=2000,per.jump=0.7){
    p<-ncol(E)
    n.pc<-ncol(U)

    # PC
    E.pc<-E%*%U

    # % of variance
    var.total<-sum(diag(cov(E)))
    var.per=var.per.ind<-rep(NA,n.pc)
    for(j in 1:n.pc)
    {
      var.per.ind[j]<-var(E.pc[,j])/var.total
      var.per[j]<-sum(diag(cov(matrix(E.pc[,1:j],ncol=j))))/var.total
    }

    lambda.est<-rep(NA,n.pc)
    V<-matrix(NA,p,n.pc)
    var.per.new<-rep(NA,n.pc)
    # first PC lambda choice
    if(is.null(D))
    {
      out.tmp<-fusedlasso1d(y=E.pc[,1],X=E,gamma=gamma,eps=eps,maxsteps=maxsteps)
    }else
    {
      #out.tmp<-fusedlasso(y=E.pc[,1],X=E,D=D,gamma=gamma,eps=eps,maxsteps=maxsteps)
      out.tmp<-genlasso(y=E.pc[,1],X=E,D=D,gamma=gamma,eps=eps,maxsteps=maxsteps)
    }
    var.per.tmp<-rep(NA,length(out.tmp$lambda))
    for(k in 1:length(out.tmp$lambda))
    {
      var.per.tmp[k]<-var(out.tmp$fit[,k])/var.total
    }
    var.per.diff.tmp<-var.per.tmp[2:length(var.per.tmp)]-var.per.tmp[1:(length(var.per.tmp)-1)]
    lambda.idx.tmp<-max(which(var.per.diff.tmp>quantile(var.per.diff.tmp,probs=per.jump)))+1
    #Compares the drop in variance explained for higher values of lambda.
    #Chooses the lowest lambda for which the percentage of variance explained
    #drops by at least the per.jump percentile of all drops
    #Since bigger drops occur for higher lambdas, a high value of per.jump
    #will favor higher lambdas and more shrinkage.

    # lambda.idx.tmp<-min(which(abs(var.per.tmp-var.per.ind[1])<var.diff))
    lambda.est[1]<-out.tmp$lambda[lambda.idx.tmp]
    V[,1]<-out.tmp$beta[,lambda.idx.tmp]
    var.per.new[1]<-var.per.tmp[lambda.idx.tmp]

    for(j in 2:n.pc)
    {
      Etmp<-deCor(cbind(E%*%V[,1:(j-1)],E.pc[,j]))

      if(is.null(D))
      {
        out.tmp<-fusedlasso1d(y=Etmp[,j],X=E,gamma=gamma,eps=eps,maxsteps=maxsteps)
      }else
      {
        #out.tmp<-fusedlasso(y=Etmp[,j],X=E,D=D,gamma=gamma,eps=eps,maxsteps=maxsteps)
        out.tmp<-genlasso(y=Etmp[,j],X=E,D=D,gamma=gamma,eps=eps,maxsteps=maxsteps)
      }

      var.per.tmp=var.per.tol.tmp<-rep(NA,length(out.tmp$lambda))
      for(k in 1:length(out.tmp$lambda))
      {
        var.per.tmp[k]<-var(out.tmp$fit[,k])/var.total
        dtmp<-deCor(cbind(E%*%V[,1:(j-1)],out.tmp$fit[,k]))
        var.per.tol.tmp[k]<-sum(diag(cov(dtmp)))/var.total
      }
      var.per.diff.tmp<-var.per.tol.tmp[2:length(var.per.tol.tmp)]-var.per.tol.tmp[1:(length(var.per.tol.tmp)-1)]
      lambda.idx.tmp<-max(which(var.per.diff.tmp>quantile(var.per.diff.tmp,probs=per.jump)))+1
      lambda.est[j]<-out.tmp$lambda[lambda.idx.tmp]
      V[,j]<-out.tmp$beta[,lambda.idx.tmp]
      var.per.new[j]<-var.per.tol.tmp[lambda.idx.tmp]
    }

    re<-list(lambda=lambda.est,V=V,var.per=var.per.new)
    return(re)
  }

mediate_multiple <- #Fixed!
  function(X, M, Y, sims = 1000, boot.ci.type = c("bca", "perc"),
           conf.level = 0.95)
  {
    n <- nrow(M)
    p <- ncol(M)

    if (is.null(colnames(M))){
      colnames(M) <- paste0("m", 1:p)
    }

    coln <- c("estimate", "se", "pvalue", "lbound", "ubound")

    #Empty matrices to fill with results
    alpha <- matrix(NA, p, 5)
    colnames(alpha) <- coln
    rownames(alpha) <- colnames(M)
    beta <- alpha
    beta_a <- alpha
    IE <- alpha

    IE.sims <- matrix(NA, sims, p)
    TE.sims <- matrix(NA, sims, p)
    DE.sims <- rep(NA, sims)
    TE.est <- rep(NA, sims)
    TE.pv <- rep(NA, sims)

    #Evaluate each PC as a mediator on its own
    for(j in 1:p){
      dat.tmp <- data.frame(X = X, M = M[, j], Y = Y)

      fit.m <- lm(M ~ X, data = dat.tmp)
      fit.y <- lm(Y ~ X + M, data = dat.tmp)

      alpha[j, 1:5] <-
        c(coef(fit.m)[2],
          summary(fit.m)$coefficients[2,2],
          summary(fit.m)$coefficients[2,4],
          confint(fit.m,level=conf.level)[2,],
        )

      beta[j, 1:5] <-
        c(coef(fit.y)[3],
          summary(fit.y)$coefficients[3,2],
          summary(fit.y)$coefficients[3,4],
          confint(fit.y,level=conf.level)[3,]
        )

      set.seed(123) #using the same seed for each makes the sims more comparable
      re.med <- mediate(fit.m, fit.y, treat = "X", mediator = "M", sims = sims,
                        boot = T, boot.ci.type = boot.ci.type[1],
                        conf.level = conf.level)

      IE[j, 1:5] <- c(re.med$d1, sd(re.med$d1.sims), re.med$d1.p, re.med$d1.ci)
      beta_a[j, 1:5] <- c(re.med$z1, sd(re.med$z1.sims), re.med$z1.p, re.med$z1.ci)

      TE.est[j] <- re.med$tau.coef #estimated total effect
      TE.pv[j] <- re.med$tau.p
      IE.sims[, j] <- re.med$d1.sims
      TE.sims[, j] <- re.med$d1.sims + re.med$z1.sims
    }

    #Organize sims
    TE.sims <- rowMeans(TE.sims)
    IE.total.sims <- rowSums(IE.sims)
    DE.sims <- TE.sims - IE.total.sims

    #Total effects
    TE <- matrix(NA,1,5)
    colnames(TE) <- coln
    TE[1,1] <- mean(TE.est) #should all be the same anyways, because seed is set
    TE[1,2] <- sd(TE.sims)
    TE[1,3] <- mean(TE.pv)

    #Global indirect effect
    IE.total <- TE
    IE.total[1,1] <- sum(IE[,1])
    IE.total[1,2] <- sd(IE.total.sims)
    IE.total[1,3] <- 2 * pnorm(abs(IE.total[1,1]/IE.total[1,2]), lower.tail=FALSE)

    #Direct effect
    DE <- TE
    DE[1,1] <- TE[1,1] - IE.total[1,1]
    DE[1,2] <- sd(DE.sims)
    DE[1,3] <- 2 * pnorm(abs(DE[1,1] / DE[1,2]), lower.tail=FALSE)

    if(boot.ci.type[1] == "bca"){

      TE[1, c(4,5)] <-
        BC.CI(TE.sims, sims = sims, conf.level = conf.level)

      IE.total[1, c(4,5)] <-
        BC.CI(IE.total.sims, sims = sims, conf.level = conf.level)

      DE[1, c(4,5)] <-
        BC.CI(DE.sims, sims = sims, conf.level = conf.level)

    } else{

      percentiles <- c((1 - conf.level) /2, 1 - (1 - conf.level) / 2)
      TE[1, c(4,5)] <- quantile(TE.sims, probs = percentiles)
      IE.total[1, c(4,5)] <- quantile(IE.total.sims, probs = percentiles)
      DE[1, c(4,5)] <- quantile(DE.sims, probs = percentiles)
    }

    re <- list(IE = IE, DE = DE, alpha = alpha, beta = beta, beta_a = beta_a,
               IE.total = IE.total, TE = TE)
    return(re)
  }

BC.CI <-
  function(theta, sims, conf.level = 0.95)
  {
    low <- (1 - conf.level)/2
    high <- 1 - low
    z.inv <- length(theta[theta < mean(theta)])/sims
    z <- qnorm(z.inv)
    U <- (sims - 1) * (mean(theta) - theta)
    top <- sum(U^3)
    under <- (1/6) * (sum(U^2))^{3/2}
    a <- top/under
    lower.inv <- pnorm(z + (z + qnorm(low))/(1 - a * (z + qnorm(low))))
    lower2 <- lower <- quantile(theta, lower.inv)
    upper.inv <- pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
    upper2 <- upper <- quantile(theta, upper.inv)
    return(c(lower, upper))
  }

deCor <-
  function(X)
  {
    n<-nrow(X)
    q<-ncol(X)

    if(q>1)
    {
      Xnew<-matrix(NA,n,q)
      Xnew[,1]<-X[,1]
      for(j in 2:q)
      {
        xtmp<-X[,1:(j-1)]
        fit<-lm(X[,j]~xtmp)
        Xnew[,j]<-fit$residuals+fit$coefficients[1]
      }
    }else
    {
      Xnew<-X
    }
    return(Xnew)
  }

decorrelate_MX <-
  function(M.tilde,X)
  {
    n<-nrow(M.tilde)
    q<-ncol(M.tilde)

    if(q>1)
    {
      Mnew<-matrix(NA,n,q)
      Mnew[,1]<-M.tilde[,1]
      for(j in 2:q)
      {
        xtmp<-M.tilde[,1:(j-1)]
        fit<-lm(M.tilde[,j]~X+xtmp)
        Mnew[,j]<-fit$residuals+cbind(rep(1,n),X)%*%fit$coefficients[c(1,2)]
      }
    }else
    {
      Mnew<-M.tilde
    }
    return(Mnew)
  }


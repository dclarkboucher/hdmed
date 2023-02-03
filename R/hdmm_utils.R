

PDM_1 =function(x,y,m,imax, tol, theta=theta, w1=w1, interval = interval, step = step)  {

  I=diag(length(w1))

  f=function(lambda1,psi,phi){
    crossprod(solve(lambda1*I + psi) %*% phi) - 1
  }

  #negative log
  g11=function(x,y,m,theta,w1,s1,s2){
    crossprod(y-theta[3]-x*theta[5]-m%*%w1*theta[4]) / (1/s1^2)
  }

  g12=function(x,y,m,theta,w1,s1,s2){
    crossprod(m%*%w1-theta[1]-x*theta[2]) / (1/s2^2)
  }


  count=0

  for (i in 1:imax) {

    mw1=m%*%w1
    J=rep(1,nrow(m))
    X1=cbind(J,x,mw1)
    X2=cbind(J,x)

    s1=sqrt( crossprod( y-X1%*%solve( crossprod(X1) )%*% crossprod(X1, y) )  / ( nrow(m)-3 ) )

    s2=sqrt(  crossprod(mw1-X2%*%solve(  crossprod(X2) )%*% crossprod(X2, mw1)) / (nrow(m)-2) )

    psi=as.numeric(1/s1^2)*crossprod(m)*(theta[4])^2
    +as.numeric(1/s2^2)*t(m)%*%m

    phi=as.numeric(1/s1^2)*crossprod(m, y-theta[3]-x*theta[5])*theta[4]
    +as.numeric(1/s2^2)*crossprod(m, theta[1] + x*theta[2])

    grid = round(2*interval/step)


    temp=numeric( step + 1 )


    for (i in 1 : (step + 1) ) {
      temp[i]=f(- (interval + grid ) + grid*i, psi = psi, phi = phi)
    }


    for (i in 1: step ){
      if ((temp[i]*temp[i+1])<0){
        l= - (interval + grid ) + grid*i
        h= - (interval + grid ) + grid*(i+1)
      }
    }



    lambda = uniroot(f, c(l, h), phi = phi, psi = psi)$root


    w1_new=solve(lambda*I+psi)%*%phi

    theta_1=optim(theta,g11,x=x,y=y,m=m,s1=s1,s2=s2,w1=w1_new, method="BFGS")$par
    theta_2=optim(theta,g12,x=x,y=y,m=m,s1=s1,s2=s2,w1=w1_new, method="BFGS")$par
    theta_new=c(theta_2[1],theta_2[2],theta_1[3],theta_1[4],theta_1[5])

    #L = vector(length=imax)
    L = - crossprod(y-theta_new[3]-x*theta_new[5]-m%*%w1*theta_new[4]) / (2*s1^2)
    - crossprod(m%*%w1-theta_new[1]-x*theta_new[2]) / (2*s2^2)

    #     print(L[i])
    #

    if (abs( crossprod(theta_new-theta) )<tol && abs( crossprod(w1_new - w1) )<tol) break

    else {theta=theta_new
    w1=w1_new}

    count=count + 1
  }


  out1=list(w1,theta,lambda)
  names(out1) <- c("w1","theta","lambda")
  return(out1)
}

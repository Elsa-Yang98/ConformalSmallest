#' Efficiency first conformal prediction for ridge regression
#'
#' @param X A N*d training matrix
#' @param Y A N*1 training vector
#' @param X0 A N0*d testing vector
#' @param lambda a sequence of penalty parameters for ridge regression
#' @param alpha miscoverage level
#' @return upper and lower prediction intervals for X0.
#' @examples
#' df=3
#' d = 5
#' n=50   #number of training samples
#' n0=10  #number of prediction points
#' rho=0.5
#' Sigma=matrix(rho,d,d)
#' diag(Sigma)=rep(1,d)
#' beta=rep(1:5,d/5)
#' X0=mvtnorm::rmvt(n0,Sigma,df)
#' X=mvtnorm::rmvt(n,Sigma,df)	#multivariate t distribution
#' eps=rt(n,df)*(1+sqrt(X[,1]^2+X[,2]^2))
#' Y=X%*%beta+eps
#' out.efcp=efcp.fun(X,Y,X0)
#' out.efcp$up
#' out.efcp$lo
#' @export
efcp_ridge=function(X,Y,X0,lambda=seq(0,100,length=100),alpha=0.1){
  n=dim(X)[1]
  d=dim(X)[2]
  n0=dim(X0)[1]
  D1=1:floor(n/2)
  D2=(floor(n/2)+1):n
  X1=X[D1,]
  X2=X[D2,]
  Y1=Y[D1]
  Y2=Y[D2]
  n2=length(Y2)

  out.gnet = glmnet::glmnet(X1,Y1,alpha=0,lambda=lambda)
  lambda = out.gnet$lambda
  RR2=matrix(0,n2,length(lambda))
  YY2=matrix(rep(Y2,length(lambda)),nrow=length(Y2))
  RR2=abs(YY2-cbind(rep(1,n2),X2)%*%coef(out.gnet))
  q2=apply(RR2,2,quantile,probs=(1-alpha)*(1+1/n2))
  index=which(q2==min(q2))[1]
  #chosen.lambda=lambda[index]

  R2=abs(Y2-cbind(rep(1,n2),X2)%*%coef(out.gnet)[,index])
  q3=quantile(R2,probs=(1-alpha)*(1+1/n2))
  up=cbind(rep(1,n0),X0)%*%coef(out.gnet)[,index]+q3
  lo=cbind(rep(1,n0),X0)%*%coef(out.gnet)[,index]-q3

  return(list(up=up,lo=lo,lambda=lambda[index]))
}


#' Validity first conformal prediction for ridge regression
#'
#' @param X A N*d training matrix
#' @param Y A N*1 training vector
#' @param X0 A N0*d testing vector
#' @param lambda a sequence of penalty parameters for ridge regression
#' @param alpha miscoverage level
#' @return upper and lower prediction intervals for X0.
#' @examples
#' df=3
#' d = 5
#' n=50   #number of training samples
#' n0=10  #number of prediction points
#' rho=0.5
#' Sigma=matrix(rho,d,d)
#' diag(Sigma)=rep(1,d)
#' beta=rep(1:5,d/5)
#' X0=mvtnorm::rmvt(n0,Sigma,df)
#' X=mvtnorm::rmvt(n,Sigma,df)	#multivariate t distribution
#' eps=rt(n,df)*(1+sqrt(X[,1]^2+X[,2]^2))
#' Y=X%*%beta+eps
#' out.vfcp=vfcp.fun(X,Y,X0)
#' out.vfcp$up
#' out.vfcp$lo
#' @export
vfcp_ridge=function(X,Y,X0,lambda=seq(0,100,length=100),alpha=0.1){
  n=dim(X)[1]
  d=dim(X)[2]
  n0=dim(X0)[1]
  D1=1:floor(n/3)
  D2=(floor(n/3)+1):floor(2*n/3)
  D3=(floor(2*n/3)+1):n
  X1=X[D1,]
  X2=X[D2,]
  X3=X[D3,]
  Y1=Y[D1]
  Y2=Y[D2]
  Y3=Y[D3]
  n2=length(Y2)
  n3=length(Y3)

  out.gnet = glmnet::glmnet(X1,Y1,alpha=0,lambda=lambda)
  lambda = out.gnet$lambda
  RR2=matrix(0,n2,length(lambda))
  YY2=matrix(rep(Y2,length(lambda)),nrow=length(Y2))
  RR2=abs(YY2-cbind(rep(1,n2),X2)%*%coef(out.gnet))
  q2=apply(RR2,2,quantile,probs=(1-alpha)*(1+1/n2))
  index=which(q2==min(q2))[1]
  #chosen.lambda=lambda[index]

  R3=abs(Y3-cbind(rep(1,n3),X3)%*%coef(out.gnet)[,index])
  q3=quantile(R3,probs=(1-alpha)*(1+1/n3))
  up=cbind(rep(1,n0),X0)%*%coef(out.gnet)[,index]+q3
  lo=cbind(rep(1,n0),X0)%*%coef(out.gnet)[,index]-q3

  return(list(up=up,lo=lo,lambda=lambda[index]))
}


#' Cross validation conformal prediction for ridge regression
#'
#' @param X A N*d training matrix
#' @param Y A N*1 training vector
#' @param X0 A N0*d testing vector
#' @param lambda a sequence of penalty parameters for ridge regression
#' @param alpha miscoverage level
#' @param nfolds number of folds
#' @return upper and lower prediction intervals for X0
#' @export
cv.fun=function(X,Y,X0,lambda=seq(0,100,length=100),nfolds=10,alpha=0.1){
  n=nrow(X)
  n0=nrow(X0)
  D1_2split=1:floor(n/2)
  D2_2split=(floor(n/2)+1):n
  X1_2split=X[D1_2split,]
  X2_2split=X[D2_2split,]
  Y1_2split=Y[D1_2split]
  Y2_2split=Y[D2_2split]
  n22=length(Y2_2split)

  out.gnet= glmnet::cv.glmnet(X1_2split,Y1_2split,alpha=0,nfolds=nfolds,lambda=lambda)
  opt_lambda=out.gnet$lambda.min
  cv.coef=coef(out.gnet,s=opt_lambda)
  R2.cv=abs(Y2_2split-cbind(rep(1,n22),X2_2split)%*%cv.coef)
  q2.cv=quantile(R2.cv,probs=(1-alpha)*(1+1/n22))

  up=as.vector(cbind(rep(1,n0),X0)%*%cv.coef)+q2.cv
  lo=as.vector(cbind(rep(1,n0),X0)%*%cv.coef)-q2.cv
  return(list(up=up,lo=lo,lambda=opt_lambda))

}

#' Conformal prediction for ridge regression, tuning parameter by minimizing the mean of the residuals
#'
#' @param X A N*d training matrix
#' @param Y A N*1 training vector
#' @param X0 A N0*d testing vector
#' @param lambda a sequence of penalty parameters for ridge regression
#' @param alpha miscoverage level
#' @return upper and lower prediction intervals for X0
#' @export
star.fun=function(X,Y,X0,lambda=seq(0,100,length=100),alpha=0.1){
  n=dim(X)[1]
  d=dim(X)[2]
  n0=dim(X0)[1]
  D1=1:floor(n/3)
  D2=(floor(n/3)+1):floor(2*n/3)
  D3=(floor(2*n/3)+1):n
  X1=X[D1,]
  X2=X[D2,]
  X3=X[D3,]
  Y1=Y[D1]
  Y2=Y[D2]
  Y3=Y[D3]
  n2=length(Y2)
  n3=length(Y3)


  out.gnet = glmnet::glmnet(X1,Y1,alpha=0,lambda=lambda)
  lambda = out.gnet$lambda
  RR2=matrix(0,n2,length(lambda))
  YY2=matrix(rep(Y2,length(lambda)),nrow=length(Y2))
  RR2=abs(YY2-cbind(rep(1,n2),X2)%*%coef(out.gnet))
  R2=apply(RR2,2,mean)
  index.star=which(R2==min(R2))[1]
  R3.star=abs(Y3-cbind(rep(1,n3),X3)%*%coef(out.gnet)[,index.star])
  q3.star=quantile(R3.star,probs=(1-alpha)*(1+1/n3))

  up=cbind(rep(1,n0),X0)%*%coef(out.gnet)[,index.star]+q3.star
  lo=cbind(rep(1,n0),X0)%*%coef(out.gnet)[,index.star]-q3.star
  return(list(up=up,lo=lo,lambda=lambda[index.star]))

}

#' Internal function used for ginverse.fun
#' @param intercept default is TRUE
#' @param lambda a vector
ginverselm.funs = function (intercept = TRUE, lambda = 0)
{

  m = length(lambda)
  for (j in 1:m)
    train.fun = function(x, y, out = NULL) {
      n = nrow(x)
      p = ncol(x)
      v = rep(1, p)
      if (!is.null(out)) {
        chol.R = out$chol.R
      }
      else {
        chol.R = vector(mode = "list", length = m)
        for (j in 1:m) {
          chol.R[[j]] = crossprod(x)
        }
      }
      beta = matrix(0, p, m)
      for (j in 1:m) {
        beta[, j] = MASS::ginv(chol.R[[j]])%*%t(x) %*% y
      }
      return(list(beta = beta, chol.R = chol.R))
    }
  predict.fun = function(out, newx) {

    return(newx %*% out$beta)
  }
  special.fun = function(x, y, out) {
    n = nrow(x)
    p = ncol(x)

    res = y - x %*% out$beta
    for (j in 1:m) {
      s = diag(x %*% MASS::ginv(out$chol.R[[j]])%*% t(x))
      res[, j] = res[, j]/(1 - s)
    }
    return(res)
  }
  active.fun = function(out) {
    p =  nrow(out$beta)
    m = ncol(out$beta)
    act.list = vector(length = m, mode = "list")
    for (i in 1:m) act.list[[i]] = 1:p
    return(act.list)
  }
  return(list(train.fun = train.fun, predict.fun = predict.fun,
              special.fun = special.fun, active.fun = active.fun))
}

#' Internal function used for ginverse.fun
#'
my.ginverselm.funs = ginverselm.funs(lambda=0)

#' Conformal prediction for linear regression
#'
#' @param x A N*d training matrix
#' @param y A N*1 training vector
#' @param x0 A N0*d testing vector
#' @param alpha miscoverage level
#' @return upper and lower prediction intervals for X0
#' @export
ginverse.fun = function(x, y, x0,alpha=0.1) {
  n = nrow(x); n0 = nrow(x0)
  out = my.ginverselm.funs$train(x,y)
  fit = matrix(my.ginverselm.funs$predict(out,x),nrow=n)
  pred = matrix(my.ginverselm.funs$predict(out,x0),nrow=n0)
  m = ncol(pred)

  x1 = x0
  #q = qt(1-alpha/2, n-d)
  q = qnorm(1-alpha/2)
  lo = up = matrix(0,n0,m)

  for (j in 1:m) {
    #sig.hat = sqrt(sum((y - fit[,j])^2)/(n-ncol(x1)))
    sig.hat = sqrt(sum((y - fit[,j])^2)/(n))
    g = diag(x1 %*% MASS::ginv(out$chol.R[[j]])%*%t(x1))
    lo[,j] = pred[,j] - sqrt(1+g)*sig.hat*q
    up[,j] = pred[,j] + sqrt(1+g)*sig.hat*q
  }

  # Return proper outputs in proper formatting
  return(list(pred=pred,lo=lo,up=up,fit=fit))
}

#' Conformal prediction for linear regression
#'
#' @param X A N*d training matrix
#' @param Y A N*1 training vector
#' @param X0 A N0*d testing vector
#' @param alpha miscoverage level
#' @return upper and lower prediction intervals for X0
#' @export
naive.fun=function(X,Y,X0,alpha=0.1){
  x=X
  y=Y
  x0=X0
  n = nrow(x); n0 = nrow(x0)
  beta=MASS::ginv(t(x)%*%x)%*%t(x)%*%y
  R=abs(y-x%*%beta)
  q=quantile(R,probs=(1-alpha)*(1+1/n))

  up=X0%*%beta+q
  lo=X0%*%beta-q
  return(list(up=up,lo=lo))
}






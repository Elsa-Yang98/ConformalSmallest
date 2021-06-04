context("testing for efcp and vfcp")

test_that("quantile for efcp is valid", {
  df=3
  l=100    #number of dimensions
  l.lambda=100
  lambda_seq=seq(0,200,l=l.lambda)
  d = 5
  alpha=0.1
  n=50   #number of training samples
  n0=10  #number of prediction points

  rho=0.5

  Sigma=matrix(rho,d,d)
  diag(Sigma)=rep(1,d)
  beta=rep(1:5,d/5)

  X0=mvtnorm::rmvt(n0,Sigma,df)

  X=mvtnorm::rmvt(n,Sigma,df)	#multivariate t distribution
  eps=rt(n,df)*(1+sqrt(X[,1]^2+X[,2]^2))
  Y=X%*%beta+eps

  expect_equal(length(vfcp.fun(X,Y,X0)$up),nrow(X0))
  expect_equal(length(efcp.fun(X,Y,X0)$up),nrow(X0))

})




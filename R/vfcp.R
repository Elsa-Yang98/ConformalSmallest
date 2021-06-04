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
vfcp.fun=function(X,Y,X0,lambda=seq(0,100,length=100),alpha=0.1){
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

return(list(up=up,lo=lo))
}

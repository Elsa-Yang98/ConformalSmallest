# ConformalSmallest

## This package implements two  selection  algorithms  for  conformal  prediction  regions  to  obtain  the smallest prediction set in practice; these are called “efficiency first” and “validity first” conformal prediction algorithms, EFCP and VFCP for short. For details please refer to [our paper].(https://arxiv.org/abs/2104.13871).

## Example 1: Tuning free ridge regression with conformal prediction
```
library(glmnet)
library(MASS)
library(mvtnorm)
source("ginverse.fun.R")
source("functions.R")
name=paste("linear_fm_t3",sep="")


df <- 3  #degrees of freedom
l <- 60    #number of dimensions 
l.lambda <- 100
lambda_seq <- seq(0,200,l=l.lambda)
dim <- round(seq(5,300,l=l))
alpha <- 0.1
n <- 200   #number of training samples
n0 <- 100  #number of prediction points
nrep <- 100 #number of independent trials
rho <- 0.5

cov.efcp <- len.efcp <- matrix(0,nrep,l)
cov.vfcp <- len.vfcp <- matrix(0,nrep,l)
cov.naive <- len.naive <- matrix(0,nrep,l)
cov.param <- len.param <- matrix(0,nrep,l)
cov.star <- len.star <- matrix(0,nrep,l)
cov.cv10 <- len.cv10 <- matrix(0,nrep,l)
cov.cv5 <- len.cv5 <- matrix(0,nrep,l)
cov.cvloo <- len.cvloo <- matrix(0,nrep,l)


out.efcp.up <- out.efcp.lo <- matrix(0,n0,l)
out.vfcp.up <- out.vfcp.lo <- matrix(0,n0,l)
out.naive.up <- out.naive.lo <- matrix(0,n0,l)
out.param.up <- out.param.lo <- matrix(0,n0,l)
out.star.up <- out.star.lo <- matrix(0,n0,l)
out.cv10.up <- out.cv10.lo <- matrix(0,n0,l)
out.cv5.up <- out.cv5.lo <- matrix(0,n0,l)
out.cvloo.up <- out.cvloo.lo <- matrix(0,n0,l)


for(i in 1:nrep){
  cat(i,"\n")
  for (r in 1:l){
    d <- dim[r]
    set.seed(i)
    
    Sigma <- matrix(rho,d,d)
    diag(Sigma) <- rep(1,d)
    X <- rmvt(n,Sigma,df)	#multivariate t distribution
    beta <- rep(1:5,d/5)
    eps <- rt(n,df)*(1+sqrt(X[,1]^2+X[,2]^2))
    Y <- X%*%beta+eps
    
    
    X0 <- rmvt(n0,Sigma,df)
    eps0 <- rt(n0,df)*(1+sqrt(X0[,1]^2+X0[,2]^2))
    Y0 <- X0%*%beta+eps0
    
    
    out.param <- ginverse.fun(X,Y,X0,alpha=alpha)
    out.param.lo[,r] <- out.param$lo
    out.param.up[,r] <- out.param$up    
    cov.param[i,r] <- mean(out.param.lo[,r] <= Y0 & Y0 <= out.param.up[,r]) 
    len.param[i,r] <- mean(out.param.up[,r]-out.param.lo[,r]) 
    
    
    out.efcp <- efcp.fun(X,Y,X0,lambda=lambda_seq,alpha=alpha)
    out.efcp.up[,r] <- out.efcp$up
    out.efcp.lo[,r] <- out.efcp$lo    
    cov.efcp[i,r] <- mean(out.efcp.lo[,r] <= Y0 & Y0 <= out.efcp.up[,r]) 
    len.efcp[i,r] <- mean(out.efcp.up[,r]-out.efcp.lo[,r]) 
    
    
    out.vfcp <- vfcp.fun(X,Y,X0,lambda=lambda_seq,alpha=alpha)
    out.vfcp.up[,r] <- out.vfcp$up
    out.vfcp.lo[,r] <- out.vfcp$lo    
    cov.vfcp[i,r] <- mean(out.vfcp.lo[,r] <= Y0 & Y0 <= out.vfcp.up[,r]) 
    len.vfcp[i,r] <- mean(out.vfcp.up[,r]-out.vfcp.lo[,r]) 
    
    out.naive <- naive.fun(X,Y,X0,alpha=alpha)
    out.naive.up[,r] <- out.naive$up
    out.naive.lo[,r] <- out.naive$lo    
    cov.naive[i,r] <- mean(out.naive.lo[,r] <= Y0 & Y0 <= out.naive.up[,r]) 
    len.naive[i,r] <- mean(out.naive.up[,r]-out.naive.lo[,r]) 
    
    
    out.star <- star.fun(X,Y,X0,lambda=lambda_seq,alpha=alpha)
    out.star.up[,r] <- out.star$up
    out.star.lo[,r] <- out.star$lo   
    cov.star[i,r] <- mean(out.star.lo[,r] <= Y0 & Y0 <= out.star.up[,r]) 
    len.star[i,r] <- mean(out.star.up[,r] - out.star.lo[,r]) 
    
    
    out.cv5 <- cv.fun(X,Y,X0,lambda=lambda_seq,alpha=alpha,nfolds=5)
    out.cv5.up[,r] <- out.cv5$up
    out.cv5.lo[,r] <- out.cv5$lo    
    cov.cv5[i,r] <- mean(out.cv5.lo[,r] <= Y0 & Y0 <= out.cv5.up[,r]) 
    len.cv5[i,r] <- mean(out.cv5.up[,r] - out.cv5.lo[,r])
    
  }
}

df.cov <- data.frame(dim,apply(cov.param,2,mean),apply(cov.naive,2,mean),apply(cov.vfcp,2,mean),apply(cov.star,2,mean),apply(cov.cv5,2,mean), apply(cov.efcp,2,mean))

df.len <- data.frame(dim,apply(len.param,2,mean),apply(len.naive,2,mean),apply(len.vfcp,2,mean),apply(len.star,2,mean),apply(len.cv5,2,mean), apply(len.efcp,2,mean))

save(dim,cov.param, cov.naive, cov.vfcp, cov.star, cov.cv5, cov.efcp, file = "cov100_t3.RData" )
save(dim,len.param, len.naive, len.vfcp, len.star, len.cv5, len.efcp, file = "len100_t3.RData" )

```


## Example 2: Tuning free conformal quantile regression with random forest
This output the right panal of Figure 1 in [our paper].
```
df <- 3
d <- 3
l.lambda <- 100
lambda_seq <- seq(0,200,l=l.lambda)
nset <- c(50,100,500,1000,5000)
alpha <- 0.1 #miscoverage level
n0 <- 100  #number of prediction points
nrep <- 100 #number of independent trials
rho <- 0.5

evaluations <- expand.grid(1:nrep, nset, c("efficient", "valid"))
no_eval <- nrow(evaluations)
width_mat <- cov_mat <- data.frame(number = rep(0, no_eval), 
                                   rep = evaluations[,1], 
                                   nset = evaluations[,2],
                                   method = evaluations[,3])
colnames(width_mat) <- colnames(cov_mat) <- c("number", "rep", "sample size", "method")

Sigma <- matrix(rho,d,d)
diag(Sigma) <- rep(1,d) #covariance matrix for X


for(idx in 1:nrow(evaluations)){
  set.seed(evaluations[idx, 1])
  if(idx%%1 == 0){
    print(idx)
  }
  n <- evaluations[idx, 2]  
  
  X <- rmvt(n,Sigma,df)	#multivariate t distribution
  
  eps1 <- rt(n,df)*(1+sqrt(X[,1]^2+X[,2]^2))
  eps2 <- rt(n,df)*(1+sqrt(X[,1]^4+X[,2]^4))
  Y <- rpois(n,sin(X[,1])^2 + cos(X[,2])^4+0.01 )+0.03*X[,1]*eps1+25*(runif(n,0,1)<0.01)*eps2
  
  X0 <- rmvt(n0,Sigma,df)
  eps01 <- rt(n0,df)*(1+sqrt(X0[,1]^2+X0[,2]^2))
  eps02 <- rt(n0,df)*(1+sqrt(X0[,1]^4+X0[,2]^4))
  Y0 <- rpois(n0,sin(X0[,1])^2 + cos(X0[,2])^4+0.01 )+0.03*X0[,1]*eps01+25*(runif(n0,0,1)<0.01)*eps02
  
  width_mat[idx,3] <- cov_mat[idx, 3] <- n
  method <- evaluations[idx, 3]
  width_mat[idx,4] <- cov_mat[idx, 4] <- method
  width_mat[idx, 2] <- cov_mat[idx, 2] <- evaluations[idx, 1]
  
  if(method == "valid"){
    split <- c(1/2, 1/2)
  } else {
    split <- 1/2
  }  
  
  beta_grid <- seq(1e-03, 4, length = 20)*alpha
  mtry_grid <- unique(ceiling(seq(1/10, 1, length = 20)*d))
  ntree_grid <- seq(50, 400, by = 50)
  
  tmp <- try(conf_CQR_reg(X, Y, split = split, beta_grid, mtry_grid, ntree_grid, method = method, alpha = alpha))
  
  while (class(tmp)=="try-error"){
    
    tmp <- try(conf_CQR_reg(X, Y, split = split, beta_grid, mtry_grid, ntree_grid, method = method, alpha = alpha),silent=TRUE)
    
  }
  width_mat[idx, 1] <- tmp$width
  cov_mat[idx, 1] <- mean(tmp$pred_set(X0, Y0))
}



par(mfrow <- c(1,2))
width_efcp <- width_vfcp <- sd_width_efcp <- sd_width_vfcp <- NULL
#sd_efcp <- sd_vfcp <- NULL
for(n in nset){
  TMP <- width_mat[evaluations[,3] == "efficient", ]
  TMP_prime <- TMP[TMP[,3] == n,]
  
  TMP <- width_mat[evaluations[,3] == "valid", ]
  TMP_prime_vfcp <- TMP[TMP[,3] == n,]
  TMP_prime_vfcp_clean =TMP_prime_vfcp[ TMP_prime_vfcp[,1]<=10^5,1]
  
  
  width_efcp <- c(width_efcp, mean(TMP_prime[,1] / TMP_prime_vfcp[,1]))
  sd_width_efcp <- c(sd_width_efcp, sd(TMP_prime[,1]/ TMP_prime_vfcp[,1])/sqrt(nrep))
  #sd_efcp  = c(sd_efcp , sd(TMP_prime[,1])/sqrt(nrep) )
  
  width_vfcp <- c(width_vfcp, mean(TMP_prime_vfcp[,1] / TMP_prime_vfcp[,1]))
  sd_width_vfcp <- c(sd_width_vfcp, sd(TMP_prime_vfcp[,1]/ TMP_prime_vfcp[,1])/sqrt(nrep))
  #sd_vfcp  = c(sd_vfcp , sd(TMP_prime_vfcp[,1])/sqrt(nrep) )
  
}




#plot(dim, width_efcp, type = 'l', ylim = range(c(width_efcp+sd_efcp)), lwd = 2)
plot(nset, width_efcp, type = 'l', ylim =c(-10,25), lwd = 2)
lines(nset, width_efcp - sd_width_efcp, type = 'l', lty = 2, lwd = 2)
lines(nset, width_efcp + sd_width_efcp, type = 'l', lty = 2, lwd = 2)
lines(nset, width_vfcp, type = 'l', ylim = range(c(width_efcp, width_vfcp)), lwd = 2, col = "red")
lines(nset, width_vfcp - sd_width_vfcp, type = 'l', lty = 2, lwd = 2, col = "red")
lines(nset, width_vfcp + sd_width_vfcp, type = 'l', lty = 2, lwd = 2, col = "red")
abline(h = 1)

cov_efcp <- cov_vfcp <-sd_cov_efcp <- sd_cov_vfcp <- NULL
for(n in nset){
  TMP <- cov_mat[evaluations[,3] == "efficient", ]
  TMP_prime <- TMP[TMP[,3] == n,]
  cov_efcp <- c( cov_efcp, mean(TMP_prime[,1] ) )
  sd_cov_efcp <- c(sd_cov_efcp, sd(TMP_prime[,1])/sqrt(nrep))
  
  TMP <- cov_mat[evaluations[,3] == "valid", ]
  TMP_prime <- TMP[TMP[,3] == n,]
  cov_vfcp <- c(cov_vfcp, mean(TMP_prime[,1]))
  sd_cov_vfcp <- c(sd_cov_vfcp, sd(TMP_prime[,1])/sqrt(nrep))
}
plot(nset, cov_efcp, type = 'l', ylim = c(0, 1), lwd = 2)
lines(nset, cov_vfcp, type = 'l', col = "red", lwd = 2)
abline(h = 1-alpha)

save(nset,nrep,width_mat, cov_mat, evaluations, width_efcp, sd_cov_efcp, sd_width_efcp,width_vfcp, sd_cov_vfcp,sd_width_vfcp, cov_efcp, cov_vfcp, alpha, file = "pois-100-repetitions.RData" )
```

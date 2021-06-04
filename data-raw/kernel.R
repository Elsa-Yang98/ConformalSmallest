kernel=as.matrix(read.csv("data-raw/sgemm_product.csv",header=FALSE))
X_full=kernel[,1:14]
Y_full=apply(kernel[,15:18],1,mean)
save(kernel,file='data/kernel.rda', compress='xz')

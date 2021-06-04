concrete=as.matrix(read.csv("data-raw/concrete.csv",header=FALSE,sep=","))
X_full=concrete[,1:8]
Y_full=concrete[,9]
save(concrete,file='data/concrete.rda', compress='xz')

superconduct=as.matrix(read.csv("data-raw/train.csv",header=FALSE))
X_full=superconduct[,1:81]
Y_full=superconduct[,82]
save(superconduct,file='data/superconduct.rda', compress='xz')

protein=as.matrix(read.csv("data-raw/CASP.csv",header=TRUE))
X_full=protein[,2:9]
Y_full=protein[,1]
save(protein,file='data/protein.rda', compress='xz')

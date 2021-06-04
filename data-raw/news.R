news=as.matrix(read.csv("data-raw/OnlineNewsPopularity.csv",header=FALSE))
X_full=news[,1:59]
Y_full=news[,60]
save(news,file='data/news.rda', compress='xz')

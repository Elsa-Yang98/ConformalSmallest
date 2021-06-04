blog=as.matrix(read.csv("data-raw/blogData_train.csv",header=TRUE))
X_full=blog[,1:280]
Y_full=blog[,281]
save(blog,file='data/blog.rda', compress='xz')

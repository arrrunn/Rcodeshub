setwd("D:/Academics - Rutgers/2. PhD/2. Spring 2020/597 Data Wrangling/HW3")
data.train <- read.table('vowel.train.txt',sep=",",header=T)
data.test <- read.table('vowel.test.txt',sep=",",header=T)

################ User functions ###################################################

LDA <- function(x) {
  #  t <- t(x[1,])%*% solve(covar) %*% t(mu)  - diag(mu %*% solve(covar) %*% t(mu))
  t <- t(apply(x,1,function(x) {t(x)%*% solve(covar) %*% t(mu)  - 1/2*diag(mu %*% solve(covar) %*% t(mu)) + t(log(pi))}))
  class <- apply(t,1,function(x) {which(x==max(x))})
  return(class)
}


QDA <- function(x) {
  t=NULL
  #  t <- t(x[1,])%*% solve(covar) %*% t(mu)  - diag(mu %*% solve(covar) %*% t(mu))
  for (k in 1:length(cat)) {
  t1 <- apply(x,1,function(x) {t(log(pi[k])) - 1/2*log(det(covar.k[[k]]))  - 1/2*(t(x-mu[k,]) %*% solve(covar.k[[k]]) %*% (x-mu[k,])) })
  t <- cbind(t,t1)
  }
  colnames(t) <- c(1:length(cat))
  class <- apply(t,1,function(x) {which(x==max(x))})
  return(class)
}


###################################################################################
#Problem 1
data.train <- data.train[-1]
data.test  <- data.test [-1]

cat <- levels(factor(data.train$y))

#Estaimted values
mu <- (apply(data.train[,-1],2,function(x) {tapply(x,data.train$y,mean)})) #k x p matrix
pi <- as.matrix(tapply(data.train[,1],data.train$y,length)) #k x 1 vector
#table(data.train$y)

x_norm <- t(apply(data.train,1,function(x) {x[-1]-mu[x[1],]})) #n x p matrix
#data.train[1,-1]-mu[1,] Check

#LDA
covar <- 1/(dim(x_norm)[1]-length(cat)) * t(x_norm) %*% x_norm #since N's are equal, just divide covariance matrix by k

pred_LDA <- LDA(as.matrix(data.train[-1]))
ct <- table(data.train$y,pred_LDA)
ct
sum(diag(prop.table(ct)))

pred_LDA <- LDA(as.matrix(data.test[-1]))
ct <- table(data.test$y,pred_LDA)
ct
sum(diag(prop.table(ct)))


#QDA
covar.k <- NULL
for (i in 1:length(cat)) {
covar.k[i] <- list(1/(length(data.train$y==i)-1)*t(x_norm[data.train$y==i,]) %*% x_norm[data.train$y==i,])
}

pred_QDA <- QDA(as.matrix(data.train[-1]))
ct <- table(data.train$y,pred_QDA)
ct
sum(diag(prop.table(ct)))

pred_QDA <- QDA(as.matrix(data.test[-1]))
ct <- table(data.test$y,pred_QDA)
ct
sum(diag(prop.table(ct)))

###################################################################################
# Problem 2
x1 <- rnorm(25,3,1)
x2 <- rnorm(25,-4,1)

y1 <- rnorm(25,0,1)
y2 <- rnorm(25,0,1)

library(cluster)
data <- rbind(cbind(x1,x2),cbind(y1,y2))
plot(data)
k1 <- kmeans(data,2)
k2 <- kmeans(data,3)
par(mfrow=c(1,2))
plot(data,col=k1$cluster)
plot(data,col=k2$cluster)

d <- dist(data)
h1 <- hclust(d,method="single")
h2 <- hclust(d,method="complete")
h3 <- hclust(d,method="average")

par(mfrow=c(1,3))
plot(h1,main="Single Linkage")
plot(h2,main="Complete Linkage")
plot(h3,main="Average Linkage")

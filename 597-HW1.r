install.packages("DescTools")
library(DescTools)

setwd("D:/Academics - Rutgers/2. PhD/2. Spring 2020/597 Data Wrangling/HW1")
data.train <- read.table('zip.train.txt',sep=" ")
data.test <- read.table('zip.test.txt',sep=" ")

table(data.train$V1)
table(data.test$V1)

data.train$V258 <- NULL

dat1 <- data.frame(apply(data.train,2,as.numeric))
mod1 <- lm(V1 ~ .,data=dat1)
summary(mod1)

dat1.pred <- cbind(data.test[1],predict.lm(mod1,data.test[-1]))
plot(dat1.pred[,1],dat1.pred[,2])
dat1.pred[,3] <- round(dat1.pred[,2])


#k-NN
k=10

dat2 <- data.frame(apply(data.train,2,as.numeric))
Y <- dat2$V1
dat2 <- dat2[-1]

summary(dat2) #to check if normalized

dat2.dist <- as.matrix(dist(dat2, method = "euclidean"))

dat2.test <- data.test[-1]
x<-dat2.test[1,]-dat2 
y<-t(dat2.test[1,]-dat2[1,])

dat3 <- sweep(dat2,2,dat2.test[1,])

function(){
  
  sum((x-y)^2)
  
}

dat2.k <- apply(dat2.dist,1,function(x) {which(Large(x,k))})


dat2.k1 <- t(apply(dat2.dist,1,function(x) {grep(paste0(Large(x,k),collapse="|"),x)}))

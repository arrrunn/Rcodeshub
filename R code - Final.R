setwd("D:/Academics - Rutgers/2. PhD/2. Spring 2020/597 Data Wrangling/Final Exam")

set.seed(123)
mu1 = runif(1,-3,1); mu2 = runif(1,0,3); prob = runif(1,0.2,0.5); sigma1 = 1; sigma2=1
print(c(mu1,mu2,prob))
X1 <- rnorm(1000, mu1, sigma1)
X2 <- rnorm(1000, mu2, sigma2)

G <- rbinom(1000,1,prob)
x <- G*X1 + (1-G)*X2

#EM function
EM.fn <- function(){
  k=2  
  while(abs(Q[k]-Q[k-1])>1E-16) {
    #E-step
    mix1 <- pi1Hat*dnorm(x, mu1Hat, sigma1Hat)
    mix2 <- pi2Hat*dnorm(x, mu2Hat, sigma2Hat)
    
    p1 <- mix1/(mix1 + mix2)
    p2 <- mix2/(mix1 + mix2)
    
    #M-Step
    pi1 <- sum(p1) / length(x)
    pi2 <- sum(p2) / length(x)
    
    mu1Hat <- sum(p1 * x) / sum(p1)
    mu2Hat <- sum(p2 * x) / sum(p2)
    
    sigma1 <- sqrt(sum.finite(p1 * (x-mu1Hat)^2) / sum(p1))
    sigma2 <- sqrt(sum.finite(p2 * (x-mu2Hat)^2) / sum(p2))
    
    p1Hat <- pi1 
    p2Hat <- pi2
    
    k=k+1
    Q[k] <- sum(log(mix1+mix2))
  }
  
  return(c(p1Hat,p2Hat,mu1Hat,mu2Hat,sigma1Hat,sigma2Hat))
}

#Initialize
Q <- 0
km1 <- kmeans(x,2)$cluster
km <- km1
if(sum(km1==1)/length(km1)>0.5){
  km[km1==2]=1
  km[km1==1]=2
  }
pi1Hat <- sum(km==1)/length(km)
pi2Hat <- sum(km==2)/length(km)
mu1Hat <- mean(x[km==1])
mu2Hat <- mean(x[km==2])
sigma1 <- sd(x[km==1])
sigma2Hat <- sd(x[km==2])

Q[1] <- sum.finite(log(pi1)*log(dnorm(x, mu1Hat, sigma1Hat))) +
          sum.finite(log(pi2)*log(dnorm(x, mu2Hat, sigma2Hat)))
Q[2] <- sum.finite(log(pi1)*log(dnorm(x, mu1Hat, sigma1Hat))) +
          sum.finite(log(pi2)*log(dnorm(x, mu2Hat, sigma2Hat)))

#Running the EM function
EM.out <- EM.fn()

pi1Pred <- EM.out[1]
pi2Pred <- EM.out[2]
mu1Pred <- EM.out[3]
mu2Pred <- EM.out[4]
sigma1Pred <- EM.out[5]
sigma2Pred <- EM.out[6]

x1<-seq(from=range(x)[1], to=range(x)[2], length.out=1000)
hist(x, prob=T,ylim=c(0,0.5))
lines(density(x), col="green", lwd=2) 
lines(x1,pi1Pred*dnorm(x1, mu1Pred, sigma1Pred)+
        pi2Pred*dnorm(x1, mu2Pred, sigma2Pred),lty=2,col="Red") 
lines(x1,0.5*dnorm(x1, mu1Pred, sigma1Pred),lty=2,col="gray") 
lines(x1,0.5*dnorm(x1, mu2Pred, sigma2Pred),lty=2,col="gray")



################3 
#Q.2: Logistic Regression vs LDA

logit.fn <- function(d,p,prob){
  
  #Apply Logistic
  mod.logit <- glm(G[1:500]~.-1,data=Train.X,family="binomial")
  miscl.train <- table(G[1:500],(mod.logit$fitted.values>=prob)*1)
  miscl.test <- table(G[501:1500],(predict(mod.logit,newdata=Test.X)>=0.5)*1)
  if(dim(miscl.train)[2]==1){
    error.logit.train <- ( miscl.train[2,1])/sum(miscl.train)
  } else {
    error.logit.train <- (miscl.train[1,2] + miscl.train[2,1])/sum(miscl.train)
  }
  if(dim(miscl.test)[2]==1){
  error.logit.test <- (miscl.test[2,1])/sum(miscl.test)
  } else {
    error.logit.test <- (miscl.test[1,2] + miscl.test[2,1])/sum(miscl.test)
  }
  
  #Apply LDA
  if(p==1){
    mod.lda <- lda(G[1:500] ~ V1 , data = Train.X.LDA)
  } else {
      mod.lda <- lda(G[1:500]~.,data=Train.X.LDA)
      }
  miscl.train <- table(G[1:500],predict(mod.lda)$class)
  miscl.test <- table(G[501:1500],predict(mod.lda,newdata=Test.X.LDA)$class)
  
  if(dim(miscl.train)[2]==1){
    error.lda.train <- ( miscl.train[2,1])/sum(miscl.train)
  } else {
    error.lda.train <- (is.null(miscl.train[1,2]) + miscl.train[2,1])/sum(miscl.train)
  }
  if(dim(miscl.test)[2]==1){
    error.lda.test <- ( miscl.test[2,1])/sum(miscl.test)
    
  } else {
    error.lda.test <- (miscl.test[1,2] + miscl.test[2,1])/sum(miscl.test)
  }
  
  Z1 <- optim(par=rep(1,p+1),logit.ll)
  coeff.l2 <- rbind(mod.logit$coefficients,Z1$par)
  error <- c(error.logit.train,error.logit.test,error.lda.train,error.lda.test)
  
  return(list(error,coeff.l2))
}

#Calculate Log-likelihood
logit.ll <- function(par){
  beta=par
  ll <- t(G[1:500]) %*% log(1/(1+exp(-as.matrix(Train.X) %*% as.matrix(beta)))) +
    t(1-G[1:500]) %*% log(1-1/(1+exp(-as.matrix(Train.X) %*% as.matrix(beta))))
  return(-ll)
}

#Setting parameter values
library(MASS)
Prob <- c(0.2,0.5)
D <- c(0,1,2,3)
P <- c(1,2,5)
parm <- expand.grid(x = D, y = P, z = Prob)


error.mean <- matrix(NA,nrow=dim(parm)[1],ncol=11)
set.seed(123)
j=1
for(j in 1:dim(parm)[1]){
  error <- matrix(NA,ncol=4,nrow=1000)
for (i in 1:dim(error)[1]){
  d=parm[j,1];p=parm[j,2];prob=parm[j,3]
  
  mu1 <- c(d,rep(0,p-1))
  Sigma <- diag(p)
  X1 <- mvrnorm(n = 1500, mu1, Sigma, tol = 1e-6, empirical = FALSE)
  X2 <- mvrnorm(n = 1500, -mu1, Sigma, tol = 1e-6, empirical = FALSE)
  
  G <- rbinom(1500,1,prob)
  X <- G*X1 + (1-G)*X2
  X <- cbind(rep(1,1500),X)
  
  Train.X <- as.data.frame(X[1:500,])
  Test.X <- as.data.frame(X[501:1500,])
  Train.X.LDA <- data.frame(X[1:500,-1])
  Test.X.LDA <- as.data.frame(X[501:1500,-1])
  colnames(Test.X.LDA) = colnames(Train.X.LDA)
  if(p==1){colnames(Train.X.LDA)="V1";colnames(Test.X.LDA)="V1"}
  
  z <- logit.fn(d,p,prob)
  error[i,] <- z[[1]]
}
error.mean[j,]  <- unlist(c(parm[j,],"Mean Error" = apply(error,2,mean),"Std Error" = apply(error,2,function(x) sqrt(var(x)))))
colnames(error.mean) <- c("d","p","Prob","Mean.Train.Logit","Mean.Test.Logit","Mean.Train.Lda","Mean.Test.Lda",
                               "Var.Train.Logit","Var.Test.Logit","Var.Train.Lda","Var.Test.Lda")
}

boxplot(error.mean[,4:7])
write.csv(error.mean,"Prob2.csv")
##################
#Problem 3: Gaussian Process Estimation
set.seed(123)
X <- runif(7)
Y <- exp(-1.4*X)*cos(7*pi*X/2)
X.sort <- sort(X)
plot(X,Y)


X1 <- runif(1000)
Y1 <- exp(-1.4*X1)*cos(7*pi*X1/2)
plot(X1,Y1)

#h <- as.matrix(dist(X,diag=T,upper=T))
#cov.fn <- function(sigma, theta) sigma^2*exp(-theta*h^2)

############# MLE Estimator
GP.MLE <- function(par) {
  beta=par[1]
  sigma=par[2]
  theta=par[3]
  n=length(X)
  e <- Y-X*beta
  h <- as.matrix(dist(X,diag=T,upper=T))
  R <- exp(-theta*h^2)
  ll <- -n/2*log(sigma^2) - 1/2*log(det(R)) - 
    1/sigma^2 * t(e) %*% solve(R) %*% (e)
  return(-ll)
}

Z1 <- optim(par=c(1,1,1),GP.MLE)
Z1

est <-Z1$par
h <- as.matrix(dist(X,diag=T,upper=T))
r <- est[2]^2*exp(-est[3]*sapply(X1,function(x) {abs(x-X)^2}))
R <- est[2]^2*exp(-est[3]*h^2)
Y1_hat <- X1*est[1] + t(r)%*%solve(R)%*%(Y-X*est[1])

MLE.MSPE <- 1/1000 * t(Y1-Y1_hat)%*%(Y1-Y1_hat)

k=10
se <- sqrt(est[2]^2-apply(r,2,function(x) {t(x)%*%solve(R)%*%X}))
UCL <- Y1_hat+1.96*se
LCL <- Y1_hat-1.96*se

GP.Pred <- data.frame(cbind(X1,Y1,Y1_hat,LCL,UCL))
colnames(GP.Pred) <- c("X1","Y1","Y1_hat","LCL","UCL")
GP.Pred <- GP.Pred[order(X1),]

#Plotting
plot(GP.Pred$X1,GP.Pred$Y1, col='black',type="l",ylim=c(-3,3)) +
lines(GP.Pred$X1,GP.Pred$Y1_hat,lty='dashed',col='red',lwd=2) 
lines(GP.Pred$X1, GP.Pred$UCL,lty='dashed', col = 'black')+
lines(GP.Pred$X1, GP.Pred$LCL,lty='dashed', col = 'black')

##############
Cond.LogLik <- function(X,Y,XVar=1, YVar=1, rho=0) {
  
  m=4
  m1=4
  m2=4
  T11 <- sum(X[1:4]^2,na.rm=T)
  T22 <- sum(Y[1:4]^2,na.rm=T)
  T12 <- sum(X[1:4]*Y[1:4],na.rm=T)
  T11.m1 <- sum(X[5:8]^2,na.rm=T)
  T22.m2 <- sum(Y[9:12]^2,na.rm=T)
  
  #observed L
  loglik <- -(m+1/2*(m1+m2))*log(2*pi) -(1/2)*m*log(XVar*YVar*(1-rho^2)) - 
    (1/(2*XVar*YVar*(1-rho^2)))*(YVar*T11 + XVar*T22 - 2*rho*sqrt(XVar*YVar)*T12) -
    1/2*(m1*log(YVar)+m2*log(XVar)) -1/(2*XVar)*(T11.m1) -1/(2*YVar)*(T11.m1)
  
  return(loglik)
}


### Q1

#Assume Bivariate Gaussian Normal f(x|sigma_i^2,rho) ~ N(0,SIGMA)
#Xvar = Sigma11
#Yvar = Sigma22

X <- c(1,1,-1,-1,2,2,-2,-2,NA,NA,NA,NA)
Y <- c(1,-1,1,-1,NA,NA,NA,NA,2,2,-2,-2)
W <- cbind(X,Y)
W[is.na(W)] <- 0

EM <- function(X,Y,XVar=1, YVar=1, rho=0) {
  
  n = length(X)
  T11 <- sum(X^2,na.rm=T)
  T22 <- sum(Y^2,na.rm=T)
  T12 <- sum(X*Y,na.rm=T)
  
  for (i in 1:100){
    #Expectation:
    T11.k1 <- rho^2*(XVar/YVar)*sum(Y^2,na.rm=T) + 4*XVar
    T22.k1 <- rho^2*(YVar/XVar)*sum(X^2,na.rm=T) + 4*YVar
    
    #Maximization:
    XVar <- T11.k1/n
    YVar <- T22.k1/n
    rho <- rho/sqrt(XVar*YVar)
    
    loglik <- -n*log(2*pi) -(1/2)*n*log(XVar*YVar*(1-rho^2)) - 
      (1/(2*XVar*YVar*(1-rho^2)))*(YVar*T11 + XVar*T22 - 2*rho*sqrt(XVar*YVar)*T12)
    print(loglik)
  }
  return(loglik)
}

xaxis <- seq(-0.6,0.6,.1)
yaxis <- seq(2,3,.2)
plot(xaxis,sapply(xaxis,function(x) EM(X,Y,XVar=1, YVar=1, rho=x)),type='l',
     ylab='',xlab='')
z <- sapply(xaxis, 
            function(x) {sapply(yaxis, 
                                function(y) {EM(X,Y,XVar=y, YVar=y, rho=x)})})
persp(yaxis,xaxis,z,theta = 210, phi = 20)

EM(X,Y,XVar=8/3, YVar=8/3, rho=1/2)
EM(X,Y,XVar=8/3, YVar=8/3, rho=-1/2)
EM(X,Y,XVar=5/2, YVar=5/2, rho=0)

EM(X,Y,XVar=2, YVar=2, rho=-0.5)
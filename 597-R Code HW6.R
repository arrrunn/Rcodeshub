X <- runif(7)
Y <- exp(-1.4*X)*cos(7*pi*X/2)
plot(X,Y)

X1 <- runif(1000)
Y1 <- exp(-1.4*X1)*cos(7*pi*X1/2)
plot(X1,Y1)


#h <- as.matrix(dist(X,diag=T,upper=T))
#cov.fn <- function(sigma, theta) sigma^2*exp(-theta*h^2)

############# MLE Estimator
loglik.fn <- function(par) {
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

Z1 <- optim(par=c(1,1,1),loglik.fn)
Z1

est <-Z1$par
h <- as.matrix(dist(X,diag=T,upper=T))
r <- est[2]^2*exp(-est[3]*sapply(X1,function(x) {abs(x-X)^2}))
R <- est[2]^2*exp(-est[3]*h^2)
Y1_hat <- X1*est[1] + t(r)%*%solve(R)%*%(Y-X*est[1])
plot(X1,Y1,col='black',cex = 1) +
points(X1,Y1_hat,col='red',cex = .5)

MLE.MSPE <- 1/1000 * t(Y1-Y1_hat)%*%(Y1-Y1_hat)

############# REML Estimator
reml.fn <- function(par) {
  n=length(X)
  beta=par[1]
  sigma=1/(n-1)*par[2]
  theta=par[3]
  e <- Y-X*beta
  h <- as.matrix(dist(X,diag=T,upper=T))
  R <- exp(-theta*h^2)
  ll <- -n/2*log(sigma^2) - 1/2*log(det(R)) - 
    1/sigma^2 * t(e) %*% solve(R) %*% (e)
  return(-ll)
}
Z2 <- optim(par=c(1,1,1),reml.fn)
Z2

est <-Z2$par
h <- as.matrix(dist(X,diag=T,upper=T))
r <- est[2]^2*exp(-est[3]*sapply(X1,function(x) {abs(x-X)^2}))
R <- est[2]^2*exp(-est[3]*h^2)
Y1_hat <- X1*est[1] + t(r)%*%solve(R)%*%(Y-X*est[1])
plot(X1,Y1,col='black',cex = 1) +
  points(X1,Y1_hat,col='red',cex = .5)

REML.MSPE <- 1/1000 * t(Y1-Y1_hat)%*%(Y1-Y1_hat)

############# CV Estimator
cv.fn <- function(par) {
  CV_EMSPE = 0
  for (i in 1:length(X)){
  x = X[i]
  x1 = X[-i]
  y=Y[i]
  Y1 = Y[-i]
  theta=par
  n=length(X)
  h <- as.matrix(dist(x1,diag=T,upper=T))
  R <- exp(-theta*h^2)
  betaHat <- c(solve((t(x1)%*%solve(R)%*% x1))%*% t(x1)%*%solve(R)%*%Y1)
  sigma2Hat <- 1/n*(t(Y1-x1*betaHat) %*% solve(R) %*% (Y1-x1*betaHat))
  r <- c(sigma2Hat)^2*exp(-theta*{abs(x-x1)^2})
  Y_hat <- x*betaHat + t(r)%*%solve(R)%*%(Y1-x1*betaHat)
  CV_EMSPE <- CV_EMSPE + (Y_hat - y)^2
  }
  return(CV_EMSPE)
}

Z3 <- optim(par=1,cv.fn,method="Brent",lower=0,upper=100)
CV.MSPE=cv.fn(Z3$par)


#####################################################################

loglik.fn2 <- function(par) {
  theta=par
  n=length(X)
  h <- as.matrix(dist(X,diag=T,upper=T))
  R <- exp(-theta*h^2)
  betaHat <- c(solve((t(X)%*%solve(R)%*% X))%*% t(X)%*%solve(R)%*%Y)
  sigma2Hat <- 1/n*(t(Y-X*betaHat) %*% solve(R) %*% (Y-X*betaHat))
  ll <- n*log(sigma2Hat) + log(det(R))
  return(ll)
}
# and plot it just to make sure
xaxis <- seq(-500,500,.1)
plot(xaxis,sapply(xaxis,FUN=loglik.fn2),type='l',
     ylab='',xlab='')

loglik.fn2(0)

Z2 <- optim(1,loglik.fn2)
Z2
  
# and plot it just to make sure
xaxis <- seq(-500,500,.1)
input <- cbind(xaxis,c(rep(1,length(xaxis))),c(rep(1,length(xaxis))))
plot(xaxis,apply(input,1,FUN=loglik.fn),type='l',
     ylab='',xlab='')

cov.fn(1,2)
det(as.matrix(cov.fn(2,2)))
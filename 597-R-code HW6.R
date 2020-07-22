x <- c(0.01,0.48,0.71,0.95,1.19,0.01,0.48,1.44,0.71,1.96,0.01,1.44,1.96)
y <- c(127.6,124,110.8,103.9,101.5,130.1,122,92.3,113.1,83.7,128,91.4,86.2)

data <- cbind("x"=x,"y"=y)
mod.act <- lm(y ~ x,data=as.data.frame(data))
summary(mod.act)

#Method 1: resampling the joint distribution
coeff <- NULL
iter = 10000
for (i in 1:iter) {
  seed <- rbinom(13,13,0.5)
  boot.data <- data[seed,]
  mod <- lm(y ~ x,data=as.data.frame(boot.data))
  summary(mod)
  coeff <- cbind(coeff,as.matrix(mod$coefficients,ncol=1))
}
coeff <- t(coeff)
coeff <- cbind(coeff,'theta'=coeff[,2]/coeff[,1])
c("Mean"=mean(coeff[,3]),c(quantile(coeff[,3],0.025),quantile(coeff[,3],0.975)))
boxplot(coeff[,3])

#Method 2: resampling the residuals
coeff2 <- NULL
iter = 10000
for (i in 1:iter) {
  seed <- rbinom(13,13,0.5)
  boot.data2 <- cbind("y"=(mod.act$fitted.values + mod.act$residuals[seed]),data[,2])
  mod2 <- lm(y ~ x,data=as.data.frame(boot.data2))
  summary(mod2)
  coeff2 <- cbind(coeff2,as.matrix(mod2$coefficients,ncol=1))
}
coeff2 <- t(coeff2)
coeff2 <- cbind(coeff2,'theta'=coeff2[,2]/coeff2[,1])
c("Mean"=mean(coeff2[,3]),c(quantile(coeff2[,3],0.025),quantile(coeff2[,3],0.975)))
boxplot(coeff2[,3])

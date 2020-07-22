setwd('D:/Academics - Rutgers/2. PhD/2. Spring 2020/597 Data Wrangling/HW5')

dat <- read.csv('wage.txt',sep="")

par(mfrow=c(2,2))

plot(dat$age,dat$wage,,xlab="Age",ylab="Wages",main="Scatter Plot")
agelims<-range(dat$age)
age.grid<-seq(from=agelims[1], to = agelims[2])


#polynomial regression 
mod1 <- lm(wage ~ poly(age, 4),data=dat)
summary(mod1)
plot(dat$age,dat$wage,col="darkgrey",xlab="Age",ylab="Wages",main="Polynomial Regression") +
  points(age.grid,predict(mod1,newdata = list(age=age.grid)),col="red",lwd=2,type="l")

#cubic spline with knots at 25,40,60
library(splines)
mod2 <- lm(wage ~ bs(age,knots = c(25,40,60)),data = dat )
summary(mod2)
#Plotting the Regression Line to the scatterplot   
plot(dat$age,dat$wage,col="darkgrey",xlab="Age",ylab="Wages",main="Cubic Splines") +
  points(age.grid,predict(mod2,newdata = list(age=age.grid)),col="darkgreen",lwd=2,type="l") + 
  abline(v=c(25,40,60),lty=2,col="darkgreen")

#smoothing splines
mod3 <- smooth.spline(x=dat$age, y=dat$wage, cv=F)
plot(dat$age,dat$wage,col="darkgrey",xlab="Age",ylab="Wages",main="Smoothing Splines") +
  lines(mod3,lwd=2,col="blue")

#Plotting all three together
plot(dat$age,dat$wage,col="darkgrey",xlab="Age",ylab="Wages",main="Comparison of Methods") +
  abline(v=c(25,40,60),lty=2,col="darkgreen") +
  points(age.grid,predict(mod1,newdata = list(age=age.grid)),col="red",lwd=2,type="l",lty=1) + 
  points(age.grid,predict(mod2,newdata = list(age=age.grid)),col="darkgreen",lwd=2,type="l",lty=2) + 
  lines(mod3,lwd=2,col="blue",lty=3) +
  legend("topright", legend=c("Poly reg", "Cubic Spline", "Smoothing Spline"),
         col=c("red", "darkgreen", "blue"), lty=1:3, cex=0.8)
  
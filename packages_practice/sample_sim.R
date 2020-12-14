
#### McCulloh Tutorial

### uses BayesTree package

#simulate daa
n = 100
sigma= .1
f=function(x){x^3}
set.seed(14)
x = sort(2*runif(n)-1)
y = f(x) + sigma*rnorm(n)
xtest = seq(-1,1, by = .2)
plot(x,y)
points(xtest, rep(0, length(xtest)), col = 'red', pch = 16)

##install packages
#install.packages("BayesTree")
library(BayesTree)
rb = bart(x,y,xtest)
length(xtest)
dim(rb$yhat.test)

plot(x,y)
lines(xtest, xtest^3, col = 'blue')
lines(xtest, apply(rb$yhat.test, 2, mean), col = 'red')
qm = apply(rb$yhat.test, 2, quantile, probs = c(.05, .95))
lines(xtest, qm[1,], col = 'red', lty=2)
lines(xtest, qm[2,], col = 'red', lty = 2)



#### BAR package

## sample code from the Bayes Tree documentation

# install.packages("BART")
library(BART)

##simulate data (example from Friedman MARS paper)
f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma = 1.0  #y = f(x) + sigma*z , z~N(0,1)
n = 100      #number of observations
set.seed(99)
x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
Ey = f(x)
y=Ey+sigma*rnorm(n)
lmFit = lm(y~.,data.frame(x,y)) #compare lm fit to BART later
##run BART
set.seed(99)
bartFit = bart(x,y,ndpost=200) #default is ndpost=1000, this is to run example fast.
plot(bartFit) # plot bart fit
##compare BART fit to linear matter and truth = Ey
fitmat = cbind(y,Ey,lmFit$fitted,bartFit$yhat.train.mean)
colnames(fitmat) = c('y','Ey','lm','bart')
print(cor(fitmat))
library(BART)



### example with Boston housing data

library("MASS")
x = Boston[,c(6,13)]
y = Boston$medv
head(cbind(x,y))

par(mfrow = c(2,2))
plot(x[,1], y, xlab = "x1=rm", ylab = "y=mdev")
plot(x[,2], y, xlab = "x2=lstat", ylab = "y=mdev")
plot(x[,1], x[,2], xlab = "x1=rm", ylab = "x2=lstat")

##wbar for continuous outcomes
set.seed(99)
nd = 200
burn = 50
post = wbart(x,y, nskip = burn, ndpost = nd)

names(post)
length(post$sigma)
length(post$yhat.train.mean)

plot(poast$sigma, type = "1")
abline(v = burn, lwd = 2, col = "red")

#compare bart and linear regression
lmf = lm(y~., data.frame(x,y))
fitmar = cbind(y, post$yhat.train.mean, lmf$fitted.values)
colnames(fitmat) = c("y", "BART", "Linear")
cor(fitmat)
pairs(fitmat)

## prediction funcion with BART


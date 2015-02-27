# JAMES HYER SOLEY
# NUMERICAL METHODS
# MATH 504
# 02/22/2015
# HOMEWORK #6

# Problem #2
# Part #1

library(numDeriv)

# Create Exponential Function Function
exponential <- function(x){
  ans <- -x[1]^2*x[2]^2*exp(-((x[1]^2)/2)-((x[2]^2)/2))
  return(ans)
}

# Create Function for Gradient Norm
norm <- function(grad){
  gradnorm <- sqrt(sum(grad)^2)
  return(gradnorm)
}

# Create Gradient Function
gradient <- function(x){
  dfx.1 <- x[1]*(x[1]^2-2)*x[2]^2*exp(-(x[1]^2/2)-(x[2]^2/2))
  dfx.2 <- x[2]*(x[2]^2-2)*x[1]^2*exp(-(x[2]^2/2)-(x[1]^2/2))
  grad <- c(dfx.1,dfx.2)
  return(grad)
}

dir <- -(solve(hessian(exponential,c(.1,.1)))%*%gradient(c(.1,.1)))
dir

# Newton's Method
newton <- function(x,error){
  step <- 1
  while (norm(gradient(x)) > error){
    x.new <-  x -(solve(hessian(exponential,x))%*%gradient(x))
    x <- x.new
    step <- step + 1
  }
  return(c(x,step))
}

newton(c(.1,.1),.000001)

# Part 2
eig <- eigen(hessian(exponential,c(.1,.1)))
M <- diag(c(abs(min(eig$values))+.1,abs(min(eig$values))+.1))*diag(c(1,1))+hessian(exponential,c(.1,.1))
new.dir <- -(solve(M)%*%gradient(c(.1,.1)))
new.dir

# Newton's Method
newton <- function(x,error){
  step <- 1
  while (norm(gradient(x)) > error){
    M <- diag(c(abs(min(eig$values))+.1,abs(min(eig$values))+.1))*diag(c(1,1))+hessian(exponential,x)
    x.new <-  x -(solve(M)%*%gradient(x))
    x <- x.new
    step <- step + 1
  }
  return(c(x,step))
}

newton(c(.1,.1),.000001)

# Problem #3
# Part A
library(numDeriv)
nls <- read.table("~/Math 504/HW#6/nonlinear_least_squares.txt", header=TRUE, quote="\"")
y <- nls$y
x <- nls$x

# Create non.linear Function Function
loss <- function(d,r){
  ans <- sum((y-d*exp(-r*x))^2)
  return(ans)
}

# Create Function for Gradient Norm
norm <- function(grad){
  gradnorm <- sqrt(sum(grad)^2)
  return(gradnorm)
}

# Create Gradient Function
gradient <- function(d,r){
  df.d <- sum(-2*y*exp(-r*x)+2*d*exp(-2*r*x))
  df.r <- sum(2*x*y*d*exp(-r*x)-2*x*d^2*exp(-2*r*x))
  grad <- c(df.d,df.r)
  return(grad)
}

# Create Hessian Function
hess.loss <- function(d,r){
  drr <- sum(-2*x^2*y*d*exp(-r*x)+4*x^2*d^2*exp(-2*r*x))
  ddr <- sum(2*x*y*exp(-r*x)-4*x*d*exp(-2*r*x))
  ddd <- sum(2*exp(-2*r*x))
  hess <- matrix(c(drr,ddr,ddr,ddd),nrow=2,ncol=2)
  return(hess)
}

# Newton's Method
newton <- function(d,r,error){
  step <- 1
  eig <- eigen(hess.loss(d,r))
  print(eig)
  x <- c(d,r)
  while (norm(gradient(d,r)) > error){
    M <- diag(c(abs(min(eig$values))+.1,abs(min(eig$values))+.1))*diag(c(1,1))+hess.loss(d,r)
    print(M)
    x.new <-  x -(solve(M)%*%gradient(d,r))
    print(x.new)
    x <- x.new
    step <- step + 1
    if (step == 2){
      break
    }
  }
  return(c(d,r,step))
}

newton(1,1,.01)
debug(newton)

x <- gradient(1,1)
norm(x)

# Part B
nls.fit <- nls(y ~ d*exp(-r*x), start=list(d=2000,r=.2))
plot(x, predict(nls.fit, newx=seq(0,15,0.01)), ylim=c(0,1700))
points(x, y, col="dark gray")

summary(nls.fit)

# Question #4
rm(list=ls())
install.packages("ElemStatLearn")
library(ElemStatLearn)

head(prostate)
data(prostate)
B <- as.matrix(prostate$lpsa)
M <- as.matrix(prostate[,1:8])
head(M)

# Coordinate Descent
# scale all the predictors
predictors <- scale(prostate[,-c(9,10)])

# point magnitude function
magnitude <- function(point) {
  return(sqrt(t(point) %*% point))
}


softThreshold <- function(z, lambda) 
{
  if(z > lambda)
    return(z - lambda)
  else if(abs(z) < lambda)
    return(0)
  else if(z < -lambda)
    return(z + lambda)
}

# coordinate descent algorithm
coordinateDescent <- function(startingPoint, dependentVar, predictors, lambda, epsilon)
{
  # start out with the initial parameter values
  point <- startingPoint
  workingPoint <- point 
  nextPoint <- rep(NA, length(point))
  
  # iterate through the predictors
  while(1)
  {
    # iterate through the predictors
    for(i in 1:dim(predictors)[2])
    {
      # calculate the residuals based on all other predictors 
      r <- dependentVar - predictors[,-i] %*% workingPoint[-i]
      
      # calculate the least squares coefficient of the residuals regressed on target predictor
      fitCoef <- coefficients(lm(r ~ predictors[,i] - 1))[1]
      
      # soft thresholding step
      # store the values in the working point and the next point
      workingPoint[i] <- softThreshold(fitCoef, lambda)
      nextPoint[i] <- softThreshold(fitCoef, lambda)
    }
    
    if(magnitude(nextPoint - point) < epsilon)
      return(nextPoint) 
    
    point <- nextPoint
    workingPoint <- point 
    nextPoint <- rep(NA, length(point))
  } 
}

tail(predictors)
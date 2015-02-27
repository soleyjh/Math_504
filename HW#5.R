# JAMES HYER SOLEY
# NUMERICAL METHODS
# MATH 504
# 02/15/2015
# HOMEWORK #5

# Problem #2b
# Create function for the gradient

library(numDeriv)

# Create Banana Function Function
banana <- function(x){
  ans <- 100*(x[2]-x[1]^2)^2+(1-x[1])^2
  return(ans)
}

# Create function for gradient norm
norm <- function(grad){
  gradnorm <- sqrt(sum(grad)^2)
  return(gradnorm)
}

# Create Numerator of NM Algo
gradient <- function(x){
  dfx.1 <- 2*(200*x[1]^3-200*x[1]*x[2]+x[1]-1)
  dfx.2 <- 200*(x[2]-x[1]^2)
  grad <- c(dfx.1,dfx.2)
  return(grad)
}

# Newton's Method
newton <- function(x,error){
step <- 1
  while (norm(gradient(x)) > error){
    x.new <-  x - (solve(hessian(banana,x))%*%gradient(x))
    x <- x.new
    step <- step + 1
  }
  return(c(x,step))
}

newton(c(4,4),.000001)

# Problem #2c
# Create function for gradient norm
norm <- function(grad){
  gradnorm <- sqrt(sum(grad)^2)
  return(gradnorm)
}

# Create Numerator of NM Algo
gradient <- function(x){
  dfx.1 <- 10*(5*x^4-15*x^2+4)
  grad <- c(dfx.1)
  return(grad)
}

# Second Derivative of NW Algo for 1D
gradient.two <- function(x){
  dfx.1 <- 100*x*(2*x^2-3)
  grad <- c(dfx.1)
  return(grad)
}

# Newton's Method
newton <- function(x,error){
  step <- 1
  path <- x
  while (norm(gradient(x)) > error){
    x.new <-  x - (gradient(x)/gradient.two(x))
    path <- c(path,x.new)
    x <- x.new
    step <- step + 1
  }
  return(c(x,step))
}

newton( 3,.000001)
newton(-3,.000001)
plot(ans, ylim=c(-3,3))
lines(ans.2)

# Problem #4c
# Calculate the Eigenvalues and Vectors of A

A <- matrix(c(2,6,14,6,10,13,14,13,12), nrow=3)
A

eigen(A)
Q <- eigen(A)$vectors
Q
D <- diag(eigen(A)$values,nrow=3,ncol=3)
D

Q%*%D%*%t(Q)
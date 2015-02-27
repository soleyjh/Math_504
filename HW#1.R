# JAMES HYER SOLEY
# NUMERICAL METHODS
# MATH 611
# 01/18/2015

# Problem #3a
economic_data <- read.table("~/Math 504/HW#1/economic_data.txt", header=TRUE, quote="\"")
B <- economic_data$B
economic_data$Index <- NULL 
economic_data$B <- NULL 
economic_data$A0 <- 1

M <- as.matrix(economic_data)
M

alpha <- solve(t(M)%*%M, tol=10^-30)%*%t(M)%*%B
alpha

# Problem #3c
econ.reg<-lm(B~A1+A2+A3+A4+A5+A6, data=economic_data)
summary(econ.reg)

# Problem #3b
# Create function for the gradient
gradient <- function(x){
  df.alpha <- 2*t(M.s)%*%M.s%*%alpha - 2*t(M.s)%*%Y.s
  grad <- as.vector(df.alpha)
  return(grad)
}

# Create function for gradient norm
norm <- function(grad){
  gradnorm <- sqrt(sum(grad)^2)
  return(gradnorm)
}

# Create function for h.prime
h.prime <- function(x,s,d){
  s.new <- t(gradient(x+s*d))%*%(d)
  return(s.new)
}

# Create bisection function
bisection <- function(x){
  a <- -10
  b <- 10
  d <- gradient(x)/norm(gradient(x))
  tol <- .1
  while ((b-a)/2 > tol){
      c <- (a+b)/2
      if (h.prime(x,c,d)==0){breaK}
      if (h.prime(x,c,d)*h.prime(x,a,d) < 0){b <- c}
      else{a <- c}
    }
return(c)
}

#Create grad.descent function with all 
grad.descent <- function(starts = c(0,0,0,0,0,0,0), target = .001, step = .01){
#Initialize variables for the problem  
  x.new <- starts
  grad.new <- gradient(x.new)
  path <- matrix(0,1,7)
  gradients <- matrix(0,1,7)

while (norm(grad.new) > target) {
  if (grad.new[1] %in% gradients[,1] && grad.new[2] %in% gradients[,2]) {
    step <- step/2
 }
  path <- rbind(path, x.new)  #add new x to previous path
  x.old <- x.new 
  grad.old <- grad.new 
  gradients <- rbind(gradients,grad.old)
  dir <- grad.old/norm(grad.old)
  x.new <- x.old + dir*step
# step <- bisection(x.new)
  grad.new <- gradient(x.new)
}
return(rbind(path,x.new))
}

regression <- grad.descent()

plot(regression[,1])
plot(regression[,2])
plot(regression[,3])
plot(regression[,4])
plot(regression[,5])
plot(regression[,6])
plot(regression[,7])

tail(regression)


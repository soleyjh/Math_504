# JAMES HYER SOLEY
# NUMERICAL METHODS
# MATH 611
# 01/20/2015

# Problem #2a
# Create function for the gradient
gradient <- function(x){
  dfx.1 <- 4000*x[1]
  dfx.2 <- 2*x[2]
  grad <- c(-dfx.1,-dfx.2)
  return(grad)
}

# Create function for gradient norm
norm <- function(grad){
  gradnorm <- sqrt(sum(grad)^2)
  return(gradnorm)
}

#Create grad.descent function with all 
grad.descent <- function(starts = c(10,10), target = .01, step = .01){
  #Initialize variables for the problem  
  x.new <- starts
  grad.new <- gradient(x.new)
  path <- matrix(0,1,2)
  gradients <- matrix(0,1,2)
  
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

ans.step <- grad.descent()
tail(ans.step)
plot(ans.step[,1])
plot(ans.step[,2])

# Problem #2b
# Create function for the gradient
gradient <- function(x){
  dfx.1 <- 4000*x[1]
  dfx.2 <- 2*x[2]
  grad <- c(-dfx.1,-dfx.2)
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
  a <- 0
  b <- 20
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
grad.descent <- function(starts = c(10,10), target = .5){
  #Initialize variables for the problem  
  x.new <- starts
  grad.new <- gradient(x.new)
  path <- matrix(0,1,2)
  gradients <- matrix(0,1,2)
  step <- bisection(x.new)
  
  while (norm(grad.new) > target){
    path <- rbind(path, x.new)
    x.old <- x.new 
    grad.old <- grad.new 
    gradients <- rbind(gradients,grad.old)
    dir <- grad.old/norm(grad.old)
    x.new <- x.old + dir*step
    step <- bisection(x.new)
    grad.new <- gradient(x.new)
  }
  return(rbind(path,x.new))
}

ans.bisect <- grad.descent()
tail(ans.bisect)
plot(ans.bisect[,1])
plot(ans.bisect[,2])

# Problem #3a
# Plot All Heights on the Number Line
childrens_heights <- read.delim("~/Math 504/HW#2/childrens_heights.txt")
stripchart(childrens_heights$Boys_2[seq.int(1,26,1)], col= "red", xlim=c(0,200), xlab = "Height", main ="LDA")
stripchart(childrens_heights$Boys_9[seq.int(1,26,1)], col= "blue", add = TRUE)
stripchart(childrens_heights$Boys_18[seq.int(1,26,1)], col = "green", add = TRUE)

# Calculate Values for LDA
var.boys.2 <- var(childrens_heights$Boys_2[seq.int(1,26,1)])
var.boys.9 <- var(childrens_heights$Boys_9[seq.int(1,26,1)])
var.boys.18 <- var(childrens_heights$Boys_18[seq.int(1,26,1)])
mean.boys.2 <- mean(childrens_heights$Boys_2[seq.int(1,26,1)])
mean.boys.9 <- mean(childrens_heights$Boys_9[seq.int(1,26,1)])
mean.boys.18 <- mean(childrens_heights$Boys_18[seq.int(1,26,1)])
pooled.boys <- mean(c(var.boys.2,var.boys.9,var.boys.18))

# Create x vector with all the boys Heights
x <- c(childrens_heights$Boys_2[seq.int(1,26,1)],childrens_heights$Boys_9[seq.int(1,26,1)],childrens_heights$Boys_18[seq.int(1,26,1)])

# Create Three Class Functions
boys.2 <- function(x){
  ans <- dnorm(x, mean = mean.boys.2, sd = sqrt(pooled.boys))
  return(ans)}

boys.9 <- function(x){
  ans <- dnorm(x, mean = mean.boys.9, sd = sqrt(pooled.boys))
  return(ans)}

boys.18 <- function(x){
  ans <- dnorm(x, mean = mean.boys.18, sd = sqrt(pooled.boys))
  return(ans)}

all <-NULL
for (i in 1:length(x)){
  x.2 <- boys.2(x[i])
  x.9 <- boys.9(x[i])
  x.18 <- boys.18(x[i])
  x.new <- c(x.2,x.9,x.18)
  max <- max(x.new)
  class <- match(max,x.new)
  all <- c(all,class)
}

all

# Problem #3b
# Classify all boys and girls using LDA
# Plot All Heights on the Number Line
childrens_heights <- read.delim("~/Math 504/HW#2/childrens_heights.txt")
stripchart(childrens_heights$Boys_2[seq.int(1,26,1)], col= "red", xlim=c(0,200), xlab = "Height", main ="LDA")
stripchart(childrens_heights$Boys_9[seq.int(1,26,1)], col= "blue", add = TRUE)
stripchart(childrens_heights$Boys_18[seq.int(1,26,1)], col = "green", add = TRUE)
stripchart(childrens_heights$Girls_2[seq.int(1,32,1)], col= "orange", add = TRUE)
stripchart(childrens_heights$Girls_9[seq.int(1,32,1)], col= "yellow", add = TRUE)
stripchart(childrens_heights$Girls_18[seq.int(1,32,1)], col = "purple", add = TRUE)

# Calculate Values for LDA
var.girls.2 <- var(childrens_heights$Girls_2[seq.int(1,26,1)])
var.girls.9 <- var(childrens_heights$Girls_9[seq.int(1,26,1)])
var.girls.18 <- var(childrens_heights$Girls_18[seq.int(1,26,1)])
mean.girls.2 <- mean(childrens_heights$Girls_2[seq.int(1,32,1)])
mean.girls.9 <- mean(childrens_heights$Girls_9[seq.int(1,32,1)])
mean.girls.18 <- mean(childrens_heights$Girls_18[seq.int(1,32,1)])
pooled.girls <- mean(c(var.girls.2,var.girls.9,var.girls.18))

pooled.girls.boy <- ((pooled.girls)*32 + (pooled.boys)*26) / (26 + 32)

# Create x vector with all the boys & girls Heights
x <- c(childrens_heights$Boys_2[seq.int(1,26,1)],childrens_heights$Boys_9[seq.int(1,26,1)],childrens_heights$Boys_18[seq.int(1,26,1)],childrens_heights$Girls_2[seq.int(1,32,1)],childrens_heights$Girls_9[seq.int(1,32,1)],childrens_heights$Girls_18[seq.int(1,32,1)])

# Create Six Class Functions
boys.2 <- function(x){
  ans <- dnorm(x, mean = mean.boys.2, sd = sqrt(pooled.girls.boy))
  return(ans)}

boys.9 <- function(x){
  ans <- dnorm(x, mean = mean.boys.9, sd = sqrt(pooled.girls.boy))
  return(ans)}

boys.18 <- function(x){
  ans <- dnorm(x, mean = mean.boys.18, sd = sqrt(pooled.girls.boy))
  return(ans)}

girls.2 <- function(x){
  ans <- dnorm(x, mean = mean.girls.2, sd = sqrt(pooled.girls.boy))
  return(ans)}

girls.9 <- function(x){
  ans <- dnorm(x, mean = mean.girls.9, sd = sqrt(pooled.girls.boy))
  return(ans)}

girls.18 <- function(x){
  ans <- dnorm(x, mean = mean.girls.18, sd = sqrt(pooled.girls.boy))
  return(ans)}

all.2 <-NULL
for (i in 1:length(x)){
  b.2 <- boys.2(x[i])
  b.9 <- boys.9(x[i])
  b.18 <- boys.18(x[i])
  g.2 <- girls.2(x[i])
  g.9 <- girls.9(x[i])
  g.18 <- girls.18(x[i])
  x.new <- c(b.2,b.9,b.18,g.2,g.9,g.18)
  max <- max(x.new)
  class <- match(max,x.new)
  all.2 <- c(all.2,class)
}

all.2

# Problem #3c
# Create Six Class Functions
boys.2 <- function(x){
  ans <- dnorm(x, mean = mean.boys.2, sd = sqrt(var.boys.2))
  return(ans)}

boys.9 <- function(x){
  ans <- dnorm(x, mean = mean.boys.9, sd = sqrt(var.boys.9))
  return(ans)}

boys.18 <- function(x){
  ans <- dnorm(x, mean = mean.boys.18, sd = sqrt(var.boys.18))
  return(ans)}

girls.2 <- function(x){
  ans <- dnorm(x, mean = mean.girls.2, sd = sqrt(var.girls.2))
  return(ans)}

girls.9 <- function(x){
  ans <- dnorm(x, mean = mean.girls.9, sd = sqrt(var.girls.9))
  return(ans)}

girls.18 <- function(x){
  ans <- dnorm(x, mean = mean.girls.18, sd = sqrt(var.girls.18))
  return(ans)}

all.3 <-NULL
for (i in 1:length(x)){
  b.2 <- boys.2(x[i])
  b.9 <- boys.9(x[i])
  b.18 <- boys.18(x[i])
  g.2 <- girls.2(x[i])
  g.9 <- girls.9(x[i])
  g.18 <- girls.18(x[i])
  x.new <- c(b.2,b.9,b.18,g.2,g.9,g.18)
  max <- max(x.new)
  class <- match(max,x.new)
  all.3 <- c(all.3,class)
}

all.2
all.3

# Question #4a
MySqrt <- function(a){
  target <- .1
  s <- 1
  while (abs(s^2 - a) > target){
    s = s - ((s^2 - a)/ (2*s))}
    ans <- c(-s,s)
  return(ans)}

MySqrt(10)

# Question #4b

# Leaving blank for now

# Question #4c
uniroot(function(x) x^2 - 10, lower = 0 , upper = 10, tol = 1e-9)$root


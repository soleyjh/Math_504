# JAMES HYER SOLEY
# NUMERICAL METHODS
# MATH 504
# 03/22/2015

# PROBLEM #2
# PART A

actual <- rep(1,17)
approx <- NULL
h.path <- NULL

i <- 1
h <- 1

while (h >= 10^-16){
  h.path[i] <- h
  approx[i] <- (exp(0 + h) - exp(0 - h)) / 2*h
  h <- h/10
  i <- i + 1
}

error <- approx - actual
plot(log(h.path), error, type = 'b', xlab = "log(h)", main = "Error Plot" )

# PROBLEM #3
rm(list=ls())

## Function for the PDF of standard normal
pdf <- function(z){
  y <- (2*pi)^(-.5)*exp(-z^2/2)
  return(y)}

## Function for Riemann Sum

riemann <- function(a=a, b=b, n=n){
  dx <- abs(a-b)/n
  r <- rep(0, n)
  r[1] <- pdf(a)*dx
  
  for (i in 1:(n-1)){
    c <- a + i*dx
    r[i+1] <- pdf(c)*dx
  }
  return(sum(r))
}

## Function for Trapezoid Method

trapezoid <- function(a=a, b=b, n=n){
  dx <- abs(a-b)/n
  r <- rep(0, n+1)
  for (i in 1:n){
    c <- a + i*dx
    r[i+1] <- dx*( (pdf(c-dx) + pdf(c))/2)
  }
  return( sum( r[2:(n+1)] ) )
}

## Function for Fapprox

Fapprox <- function(a=-1.96, b= 1.96, n=100, method="reimann"){
  if(method=="reimann"){
    riemann(a=a, b=b, n=n)
  }else if(method=="trapezoid"){
    trapezoid(a=a, b=b, n=n)
  }else{    
    integrate(pdf, a, b, subdivisions=n)
  }
}

## Using a = 0 and b = 1.96, try with different values of n
options(digits = 16)
expected <- matrix(rep(.475,8), ncol=2)
results  <- matrix(rep(0,8), ncol=2)

for(i in 1:4){
  results[i,1] <- Fapprox(a = 0, b = 1.96, n = 10^i, method = "reimann")
  results[i,2] <- Fapprox(a = 0, b = 1.96, n = 10^i, method = "trapezoid")
}

error <- results - expected
error

Fapprox(a = 0, b = 1.96, n = 10, method = "UseR")
Fapprox(a = 0, b = 1.96, n = 100, method = "UseR")
Fapprox(a = 0, b = 1.96, n = 1000, method = "UseR")
Fapprox(a = 0, b = 1.96, n = 10000, method = "UseR")

# QUESTION #4
data <- read.table("~/Math_504/HW#6/nonlinear_least_squares.txt", header=TRUE, quote="\"")

#param order: d, r
fxn.4 <- list(
  f.x = function(params, data) return( sum((data$y - params[1]*exp(-params[2] * data$x) )^2) ),
  grad = function(params, data, f.x, h) { 
    d.d <- (f.x(params = c(params[1] + h, params[2]), data) - f.x(params, data)) / h
    d.r <- (f.x(params = c(params[1], params[2] + h), data) - f.x(params, data)) / h
    return(c(d.d,d.r))
  },
  hess = function(params, data, f.x) {
    return( hessian(f.x, x=params, data=data) ) 
  },
  norm = function(x) return(sqrt(sum(x)^2)) ,
  graph = function(min, data) {
    qplot(x = data$x, y = data$y) + 
      stat_function(fun = function(x) min[1]*exp(-min[2] * x), geom="line", lwd = 1, aes(colour="(0,0)")) +
      xlab("X") + ylab("Y")
  }
)

newton <- function(start, err, fxn=fxn.4, h = 10^-3, back = TRUE, grapher = FALSE, data=data, adj = .001) {
  out <- list("min" = NULL, "steps" = 0, "stepsize" = NULL, "start" = start, "path" = start, "graph" = NULL)
  x.i <- out$start
  
  while(fxn$norm(fxn$grad(x.i, data, fxn$f.x, h)) > err) {
    step <- 1
    M <- fxn$hess(x.i, data, fxn$f.x)
    
    #checking to see if any eigenvalues are negative, and if negative, shift by some lambda. 
    lambda <- min(eigen(M)$values)
    if(lambda < 0) M <- (abs(lambda) + adj)*diag(rep(1,nrow(M))) + M 
    
    #calculate direction with multiple dimensions
    dir <- -solve(M)%*%fxn$grad(x.i, data, fxn$f.x, h)
    
    #backtracking, can be turned on and off
    if(back) while(fxn$f.x((x.i + step*dir), data) > fxn$f.x(x.i, data)) step = step/2
    
    x.i <- x.i + step*dir
    out$path <- rbind(out$path, as.numeric(x.i) )
    out$steps <- out$steps + 1
    out$stepsize <- rbind(out$stepsize, step)
  }
  out$min <- as.numeric(x.i)
  
  if(grapher) out$graph <- fxn$graph(out$min, data) #graph path of results
  
  return(out)
}

newton(start = c(10,10), err = .001, fxn=fxn.4, h = 10^-4, back = TRUE, grapher = FALSE, data=data, adj = .001)
# JAMES HYER SOLEY
# NUMERICAL METHODS
# MATH 504
# 02/08/2015
# HOMEWORK #4

# Parts A & B (Testing some theories with concrete examples)
# Create Symmetric Matrix
M <- cbind(c(1,7,3),c(7,4,-5),c(3,-5,6))
M

# Create Eigenvalues of Symmetric Matrix
lambda <- eigen(M)
lambda

# Calculate the Lambda Min/Max Eigenvalues
lambda.max <- abs(lambda$values[1])
lambda.min <- abs(lambda$values[3])

# Do Test Calculations
w <- c(0,0,0)
norm.w <- sqrt(w[1]^2+w[2]^2+w[3]^2)
norm(M*w)
lambda.max*norm.w
lambda.min*norm.w

# Do Test Calculations
w <- c(5,10,2)
norm.w <- sqrt(w[1]^2+w[2]^2+w[3]^2)
norm(M*w)
lambda.max*norm.w
lambda.min*norm.w

# Do Test Calculations
w <- c(2,50,4)
norm.w <- sqrt(w[1]^2+w[2]^2+w[3]^2)
norm(M*w)
lambda.max*norm.w
lambda.min*norm.w

# Part D
A <- cbind(c(1,0,0),c(0,1,0),c(0,0,10^-21))
kappa(A)

# Compute M
M = A + t(A)

# Test for Symmetric
A == t(A)

# Compute Eigenvalues / Eigenvectors of M
eigen(M)
v.i <- eigen(M)$vectors
v.i

# Use solve() function to Break R
solve(M,v.i, tol=10^-100)
M
det(M)

# Question #3
# Part A

# Part B
economic_data <- read.table("~/Math 504/HW#1/economic_data.txt", header=TRUE, quote="\"")
B <- economic_data$B
economic_data$Index <- NULL 
economic_data$B <- NULL 
economic_data$A0 <- 1

M <- as.matrix(economic_data)
M

M.qr <- qr(M)
R <- qr.R(M.qr)
Q <- qr.Q(M.qr)

kappa(M)
kappa(M%*%t(M))
kappa(R)
kappa(Q)

alpha <- solve(R,t(Q))%*%B
alpha

# Problem #3c
econ.reg<-lm(B~A1+A2+A3+A4+A5+A6, data=economic_data)
summary(econ.reg)

# Problem #4c
# Gram-Schmidt Process

u <- cbind(as.numeric(c(3,7,5)),as.numeric(c(4,2,6)),as.numeric(c(2,3,8)))

# This is the correct answer
orthonormalization(u, basis=TRUE, norm=TRUE)

# Write my own GS function
GramSchmidt <- function(u){
  # Initialize the Matrix
  Q <- matrix(NA, dim(u)[1], dim(u)[2])
  Q[,1] <- u[,1]/norm(u[,1],"2")

  for (i in 2:ncol(u)){
    u.new <- u[,i]
  
    for (j in 1:(i-1)){
      u.new <- u.new - (Q[,j]%*%u[,i])*Q[,j]
  }  
    Q[,i] <- u.new / norm(u.new, "2")
  }
  return(Q)
}
    
GramSchmidt(u)
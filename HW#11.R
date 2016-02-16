# JAMES HYER SOLEY
# MATH 504 NUMERICAL METHODS
# 4/12/2015
# HOMEWORK #11

# Load in All Data 
d <- read.table("C:/Users/jsoley/Documents/Math_504/HW#11/images_of_three_digits.txt", sep=",")
dnba <- read.csv("C:/Users/jsoley/Documents/Math_504/HW#11/nba_ht_wt.csv",header=F,stringsAsFactor=F)
dbeer <- read.csv("C:/Users/jsoley/Documents/Math_504/HW#11/beerhall_data.csv",header=F,stringsAsFactor=F)

# PROBLEM #2
# PART A

colnames(dnba) <- c("name", "pos", "height", "weight", "age")
nba <- data.frame(dnba$height, dnba$weight)
colnames(nba) <- c("height", "weight")
plot(nba$height, nba$weight, xlab="Height", ylab="Weight", main="Height by Weight NBA")
plot(nba$height - mean(nba$height), nba$weight - mean(nba$weight), xlab="Height", ylab="Weight", main="Centered Height by Weight NBA")

# PART B
# Use Power Function From Pervious HW

norm <- function(x){
  norm <- sqrt(sum(x^2))
  return(norm)
}

power_method <- function(M, v1_old, v2_old, eps) {
  # Initialize
  N <- length(v1_old)
  v1_new <- M %*% v1_old
  v1_new <- v1_new / norm(v1_new)
  v2_new <- M %*% v2_old
  v2_new <- v2_new - as.numeric( (t(v1_new) %*% v2_new) / (t(v1_new) %*% v1_new) ) * v1_new
  v2_new <- v2_new/norm(v2_new)
  steps <- 1

  # Loop
  while(norm(v1_new - v1_old) > eps) {
    v1_old <- v1_new
    v2_old <- v2_new
    v1_new <- M %*% v1_old
    v1_new <- v1_new / norm(v1_new)
    v2_new <- M %*% v2_old
    v2_new <- v2_new <- v2_new - as.numeric( (t(v1_new) %*% v2_new) / (t(v1_new) %*% v1_new) ) * v1_new
    v2_new <- v2_new/norm(v2_new)
    steps <- steps + 1
  }

  lambda1 = t(v1_new) %*% M %*% v1_new
  lambda2 = t(v2_new) %*% M %*% v2_new
  list(v1_new = v1_new, lambda1 = lambda1, topnode1 = which.max(v1_new), 
       v2_new = v2_new, lambda2 = lambda2, topnode2 = which.max(v2_new),
       iterations = steps)
}

# Create Inputs for Power Function
eps <- 1e-10
v1 <- c(1,1)
v2 <- c(1,1)
nba.matrix <- cbind(nba$height - mean(nba$height), nba$weight - mean(nba$weight))
X <- cov(nba.matrix)

# Run Power Function
power_results <- power_method(X, v1, v2, eps)
q1 <- power_results$v1_new
q2 <- power_results$v2_new
e1 <- power_results$lambda1
e2 <- power_results$lambda2

# Check Answers against Eigen
ans <- rbind(q1,q2,e1,e2)
ans

eigen(X)

# PART C
vec <- cbind(q1,q2)
c.1 <- sapply(1:nrow(nba.matrix), function(x) t(nba.matrix[x,]) %*% vec[,1])
plot(c.1, main="Projection of C.1")
nba[which.min(c.1),]
nba[which.max(c.1),]

# PART D
vec <- cbind(q1,q2)
c.2 <- sapply(1:nrow(nba.matrix), function(x) t(nba.matrix[x,]) %*% vec[,2])
plot(c.2, main="Projection of C.2")
nba[which.min(c.2),]
nba[which.max(c.2),]

# PROBLEM #3
# PART A
# Read in Data
beer <- read.csv("C:/Users/jsoley/Documents/Math_504/HW#11/beerhall_data.csv",header=F,stringsAsFactor=F)
names(beer) <- c("county", "region", "regioncode", "crime", "hall", "school", "church")

# Create Norm Function
norm <- function(x){
  norm <- sqrt(sum(x^2))
  return(norm)
}

# Update Power Method from Part A for 3 Vars
power_method <- function(M, v1_old, v2_old, v3_old, eps) {
  # Initialize
  v1_new <- M %*% v1_old
  v1_new <- v1_new / norm(v1_new)
  v2_new <- M %*% v2_old
  v2_new <- v2_new - as.numeric( (t(v1_new) %*% v2_new) / (t(v1_new) %*% v1_new) ) * v1_new
  v2_new <- v2_new/norm(v2_new)
  v3_new <- M %*% v3_old
  v3_new <- v3_new - as.numeric( (t(v2_new) %*% v3_new) / (t(v2_new) %*% v2_new) ) * v2_new - as.numeric(( t(v1_new) %*% v3_new) / ( t(v1_new) %*% v1_new)) * v1_new
  v3_new <- v3_new/norm(v3_new)
  
  steps <- 1
  
  # Loop
  while( (norm(v1_new - v1_old) > eps) && (norm(v2_new - v2_old) > eps) && (norm(v3_new - v3_old) > eps) ) {
    v1_old <- v1_new
    v2_old <- v2_new
    v3_old <- v3_new
    v1_new <- M %*% v1_old
    v1_new <- v1_new / norm(v1_new)
    v2_new <- M %*% v2_old
    v2_new <- v2_new - as.numeric( (t(v1_new) %*% v2_new) / (t(v1_new) %*% v1_new) ) * v1_new
    v2_new <- v2_new/norm(v2_new)
    v3_new <- M %*% v3_old
    v3_new <- v3_new - as.numeric( (t(v2_new) %*% v3_new) / (t(v2_new) %*% v2_new) ) * v2_new - as.numeric(( t(v1_new) %*% v3_new) / ( t(v1_new) %*% v1_new)) * v1_new
    v3_new <- v3_new/norm(v3_new)
    
    steps <- steps + 1
  
  }
  
  lambda1 = t(v1_new) %*% M %*% v1_new
  lambda2 = t(v2_new) %*% M %*% v2_new
  lambda3 = t(v3_new) %*% M %*% v3_new
  list(v1_new = v1_new, lambda1 = lambda1, topnode1 = which.max(v1_new), 
       v2_new = v2_new, lambda2 = lambda2, topnode2 = which.max(v2_new),
       v3_new = v3_new, lambda3 = lambda3, topnode3 = which.max(v3_new),
       iterations = steps)
}

# Create Inputs for Power Function
# Used Elaine Ayo's Inputs for the Function b/c vector of 1s didn't work
eps <- 1e-5
v1 <- c(.7959,.0687,.6015)
v2 <- c(-.419,.779,.4659)
v3 <- c(.0466,.0665,-.0693)

beer.matrix <- cbind(beer$hall - mean(beer$hall), beer$school - mean(beer$school), beer$church - mean(beer$church))
X <- cov(beer.matrix)

# Run Power Function
power_results <- power_method(X, v1, v2, v3, eps)
q1 <- power_results$v1_new
q2 <- power_results$v2_new
q3 <- power_results$v3_new
e1 <- power_results$lambda1
e2 <- power_results$lambda2
e3 <- power_results$lambda3

# Check Answers against Eigen
ans <- rbind(q1,q2,q3,e1,e2,e3)
ans

eigen(X)

#Part B
vec <- cbind(q1,q2,3)
c.1 <- t(sapply(1:nrow(beer.matrix), function(x) t(beer.matrix[x,]) %*% vec[,1:2]))
data <- cbind(beer$crime, c.1)
qplot(data[,2], y=data[,3], size=data[,1], color=data[,1])


# Question #4
# Read in Data
threes <- read.table("C:/Users/jsoley/Documents/Math_504/HW#11/images_of_three_digits.txt", sep=",")

# Use Power Method From Question #2
norm <- function(x){
  norm <- sqrt(sum(x^2))
  return(norm)
}

power_method <- function(M, v1_old, v2_old, eps) {
  # Initialize
  N <- length(v1_old)
  v1_new <- M %*% v1_old
  v1_new <- v1_new / norm(v1_new)
  v2_new <- M %*% v2_old
  v2_new <- v2_new - as.numeric( (t(v1_new) %*% v2_new) / (t(v1_new) %*% v1_new) ) * v1_new
  v2_new <- v2_new/norm(v2_new)
  steps <- 1
  
  # Loop
  while(norm(v1_new - v1_old) > eps) {
    v1_old <- v1_new
    v2_old <- v2_new
    v1_new <- M %*% v1_old
    v1_new <- v1_new / norm(v1_new)
    v2_new <- M %*% v2_old
    v2_new <- v2_new <- v2_new - as.numeric( (t(v1_new) %*% v2_new) / (t(v1_new) %*% v1_new) ) * v1_new
    v2_new <- v2_new/norm(v2_new)
    steps <- steps + 1
  }
  
  lambda1 = t(v1_new) %*% M %*% v1_new
  lambda2 = t(v2_new) %*% M %*% v2_new
  list(v1_new = v1_new, lambda1 = lambda1, topnode1 = which.max(v1_new), 
       v2_new = v2_new, lambda2 = lambda2, topnode2 = which.max(v2_new),
       iterations = steps)
}

# Create Inputs for Power Function
eps <- 1e-10
v1 <- rep(1,256)
v2 <- rep(1,256)
three.cen <- sapply(1:256, function(x) threes[,x] - mean(threes[,x]))
X <- cov(three.cen)

# Run Power Function
power_results <- power_method(X, v1, v2, eps)
q1 <- power_results$v1_new
q2 <- power_results$v2_new
e1 <- power_results$lambda1
e2 <- power_results$lambda2

# Check Answers / Counts against Eigen
e1
eigen(X)$values[1]

e2
eigen(X)$values[2]

dim(q1)
dim(q2)
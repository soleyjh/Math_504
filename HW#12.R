# JAMES HYER SOLEY
# MATH 504 NUMERICAL METHODS
# 4/19/2015
# HOMEWORK #12

# QUESTION #2
# PART A

# CREATE MATRICES
set.seed(3)
A = matrix(sample(seq(100:999),4),2,2)
B = matrix(sample(c(-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9),4),2,2)
C = matrix(sample(seq(100:999)*-1,4),2,2)
Z = matrix(c(0,0,0,0),2,2)

M1 = cbind(Z,Z,A)
M2 = cbind(Z,B,Z)
M3 = cbind(C,Z,Z)
M = rbind(M1,M2,M3)
M

ans <- svd(M)
ans

ans$u%*%diag(ans$d)%*%t(ans$v)

# PART B
M1 <- ans$d[1]*ans$u[,1]%*%t(ans$v[,1]) 
M1

M2 <- ans$d[1]*ans$u[,1]%*%t(ans$v[,1]) + ans$d[2]*ans$u[,2]%*%t(ans$v[,2])
M2 

M3 <- ans$d[1]*ans$u[,1]%*%t(ans$v[,1]) + ans$d[2]*ans$u[,2]%*%t(ans$v[,2]) + ans$d[3]*ans$u[,3]%*%t(ans$v[,3])
M3

M4 <- ans$d[1]*ans$u[,1]%*%t(ans$v[,1]) + ans$d[2]*ans$u[,2]%*%t(ans$v[,2]) + ans$d[3]*ans$u[,3]%*%t(ans$v[,3]) + ans$d[4]*ans$u[,4]%*%t(ans$v[,4])
M4

M5 <- ans$d[1]*ans$u[,1]%*%t(ans$v[,1]) + ans$d[2]*ans$u[,2]%*%t(ans$v[,2]) + ans$d[3]*ans$u[,3]%*%t(ans$v[,3]) + ans$d[4]*ans$u[,4]%*%t(ans$v[,4]) + ans$d[5]*ans$u[,5]%*%t(ans$v[,5])
M5

M

# QUESTION #3
users <- read.table("~/Math_504/HW#12/user-shows.txt", quote="\"")
titles <- read.table("~/Math_504/HW#12/shows.txt", quote="\"")
ans <- svd(users)
plot(ans$d)

# Assign Using Approx
V <- ans$v
U <- ans$u

# Project Alex onto the User Space
C.3 <- t(sapply(1:nrow(users), function(x) matrix(as.numeric(users[x,]),nrow=1) %*% V[,1:25]))
alex <- matrix(as.numeric(users[500,]), nrow=1) %*% V[,1:25]
alex[1:100]

# Find Users who are like alex
tol = 3.0
pref <- data.frame(C.3)
alex.pref <- pref[abs(pref[,1] - alex[,1]) < tol
            & abs(pref[,2] - alex[,2]) < tol 
            & abs(pref[,3] - alex[,3]) < tol 
            & abs(pref[,4] - alex[,4]) < tol
            & abs(pref[,5] - alex[,5]) < tol
            & abs(pref[,6] - alex[,6]) < tol
            & abs(pref[,7] - alex[,7]) < tol
            & abs(pref[,8] - alex[,8]) < tol
            & abs(pref[,9] - alex[,9]) < tol
            & abs(pref[,10] - alex[,10]) < tol
            & abs(pref[,11] - alex[,11]) < tol
            & abs(pref[,12] - alex[,12]) < tol
            & abs(pref[,13] - alex[,13]) < tol
            & abs(pref[,14] - alex[,14]) < tol
            & abs(pref[,15] - alex[,15]) < tol
            & abs(pref[,16] - alex[,16]) < tol
            & abs(pref[,17] - alex[,17]) < tol
            & abs(pref[,18] - alex[,18]) < tol
            & abs(pref[,19] - alex[,19]) < tol
            & abs(pref[,20] - alex[,20]) < tol
            & abs(pref[,21] - alex[,21]) < tol
            & abs(pref[,22] - alex[,22]) < tol
            & abs(pref[,23] - alex[,23]) < tol
            & abs(pref[,24] - alex[,24]) < tol
            & abs(pref[,25] - alex[,25]) < tol ,]
#alex.pref
alex.pref2 <- as.numeric(rownames(alex.pref))
#alex.pref2

# Users 6414 & 8762 are Similar Lets Use those
prediction <- users[alex.pref2,1:100]
pred2 <- colSums(prediction)
tail(sort(pred2),5)

alex <- read.table("~/Math_504/HW#12/alex.txt", quote="\"")
sort(alex[,1:100])

# Question #4
# Read In Data
# PCA
senators <- read.table("~/Math_504/HW#12/senators_formatted.txt", header=TRUE, quote="\"")
votes <- read.table("~/Math_504/HW#12/votes_formatted.txt", header=TRUE, quote="\"")

# Get Eigen Vectors for Covriance of X
votes <- t(as.matrix(votes[2:101],))
votes.cen <- sapply(1:100, function(x) votes[,x] - mean(votes[,x]))
X <- cov(votes.cen)

# Project onto the First Eigenvector
vec <- cbind(eigen(X)$vectors[,1],eigen(X)$vectors[,2])
c.1 <- sapply(1:nrow(votes.cen), function(x) t(votes.cen[x,]) %*% vec[,1:2])

# Plot Values and Label
plot(c.1[1,],c.1[2,])
party <- as.vector(senators$party)
text(c.1[1,],c.1[2,],labels=party)

# SVD
senators <- read.table("~/Math_504/HW#12/senators_formatted.txt", header=TRUE, quote="\"")
votes <- read.table("~/Math_504/HW#12/votes_formatted.txt", header=TRUE, quote="\"")
votes <- t(as.matrix(votes[2:101],))

# Run SVD
ans <- svd(votes)
U <- ans$u
V <- ans$v

# Project
C <- sapply(1:nrow(votes), function(x) matrix(as.numeric(votes[x,]),nrow=1) %*% V[,1:2])
plot(C[1,],C[2,])
party <- as.vector(senators$party)
text(C[1,],C[2,],labels=party)
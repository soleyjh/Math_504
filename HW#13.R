# JAMES HYER SOLEY
# MATH 504 NUMERICAL METHODS
# 4/26/2015
# HOMEWORK #13

# PROBLEM #3
# PART A

# READ IN DATA FROM BLACKBOARD
iris <- read.csv("C:/Users/jsoley/Documents/Math_504/HW#13/iris.txt", header=FALSE)

# CENTER AND FORM COVARIANCE MATRIX
attach(iris)
iris.mat <- as.matrix(iris[1:4],)
iris.cen <- sapply(1:4, function(x) iris.mat[,x] - mean(iris.mat[,x]))
X <- cov(iris.cen)

# PROJECT ONTO TWO EIGENVECTORS
vec <- cbind(eigen(X)$vectors[,1],eigen(X)$vectors[,2])
c <- sapply(1:nrow(iris.cen), function(x) t(iris.cen[x,]) %*% vec[,1:2])

# CREATE 2D PLOT
plot(c[1,],c[2,], main="Projection Plot")
class <- iris[,5]
text(c[1,],c[2,], labels=class)
detach(iris)

# PART B
library(MASS)
lda.fit <- lda(t(c),iris[,5])

# Split Data into Training and Testing
set.seed(1)
train <- sample(1:nrow(iris),80)
test <- (-train)
iris.train <- iris[train,]
iris.test <- iris[test,]

# See How We Do USING LDA
lda.predict <- predict(lda.fit, iris.test)
table 
# JAMES HYER SOLEY
# NUMERICAL METHODS
# MATH 504 -- NUMERICAL METHODS
# 03/29/2015

# PROBLEM #2
D <- matrix(c(1,2,3,2,3,6,3,6,9),3,3)
power <- function(v,M){
  deck <- rep(0,nrow(M)) 
  int <- (M %*% v)/sqrt(sum((M%*%v)^2)) 
  deck <- rbind(deck,t(int)) 
  i <- 2 
  
  while(sqrt(sum((deck[i-1,]-deck[i,])^2)) > 10^-12){   
    old <- int   
    int <- (M%*%old)/sqrt(sum((M%*%old)^2))   
    deck <- rbind(deck,t(int))   
    i <- i+1  
  } 
  
ans <- deck[i,]  return(ans)}

eigen(D)

ans <- power(c(1,1,1),D)
ans[1]

# PROBLEM #3
# Part A

#import data
data <- read.delim("~/Dropbox/numericalmethods/ca-GrQc.txt", header=TRUE)

same <- which(data$SNode == data$ENode)

#there are 12 nodes that link to themselves

#create sparse matrix B1
#find all unique nodes, there are 5,242 unique numbers as expected
nodes <- unique(c(data$SNode, data$ENode))
nodes <- nodes[order(nodes)]
trans <- data.frame(id=1:length(nodes), num = nodes)

#transform the numbering in the data to reflect sequential numbers to make generating the matrix easier
data2<- merge(data, trans, by.x="SNode", by.y="num", all.x=TRUE)
names(data2)[3]<-"i"
data2<-merge(data2,trans,by.x="ENode",by.y="num",all.x=TRUE)
names(data2)[4]<-"j"

data2 <- data2[,3:4]

#construct B1
B1 <- matrix(0, nrow=5242, ncol=5242)

for (x in 1:28980) {
  B1[data2$i[x], data2$j[x]] <- 1
}

#check all edges accounted for, there should be 28,980
sum(B1)

#construct B2
B2 <- matrix(1, nrow=5242, ncol=5242)

#implement modified power method
norm <- function(x) sqrt(sum(x^2))

set.seed(693)
v <- runif(5242)

power <- function(B1 = B1, B2 = B2, v.old = v, tol = 10^-10, p = .15) {
  N <- length(v)
  v.new <- (1-p) * B1 %*% v.old + (p/N) * B2 %*% v.old
  v.new <- v.new / norm(v.new)
  i <- 0
  
  while( norm(v.new - v.old) > tol) {
    v.old <- v.new
    v.new <- (1-p) * B1 %*% v.old + (p/N) * B2 %*% v.old
    v.new <- v.new / norm(v.new)
    i <- i + 1
  }
  return(list(vector = v.new, 
              value = t(v.new) %*% ((1-p) * B1) %*% v.new + t(v.new) %*% ((p/N) * B2) %*% v.new, 
              topnode = which.max(v.new),
              iter = i
  ))
}

full.time <- list()
results <- list()

#now test a lot of p's 
p <- seq(0, 1, .15)
p[8] <- 1


full.time <- lapply(1:2, function(x) system.time(results[[x]] <- power(B1 = B1, B2 = B2, v.old = v, tol = 10^-10, p = x)) )

#use parallel processing
full.time <- mclapply(2:3, function(x) system.time(results[[x]] <- power(B1 = B1, B2 = B2, v.old = v, tol = 10^-10, p = x)), mc.silent = FALSE, mc.cores = getOption("cores",4), mc.cleanup = TRUE)

#now run power using sparse matrices
B1.sparse <- Matrix(0, nrow=5242, ncol=5242, sparse=TRUE)

for (x in 1:28980) {
  B1.sparse[data2$i[x], data2$j[x]] <- 1
}

#check all edges accounted for, there should be 28,980
sum(B1.sparse)

full.time.sp <- list()
results.sp <- list()

full.time.sp <- mclapply(2:3, function(x) system.time(results.sp[[x]] <- power(B1 = B1.sparse, B2 = B2, v.old = v, tol = 10^-10, p = x)), mc.silent = FALSE, mc.cores = getOption("cores",4), mc.cleanup = TRUE)

# Part B
power2 <- function(B1 = B1, B2 = B2, v1.old = u1, v2.old = u2, p = .15) {
  N <- length(v1.old)
  B1.p <- matrix(0, ncol=5242, nrow=5242)
  B1.p <- sapply(1:5242, function(x) (1-p)/sum(B1[x,]) * B1[x,])
  B2.p <- (p/N) * B2
  
  v1.new <- B1.p %*% v1.old + B2.p %*% v1.old
  v1.new <- v1.new / norm(v1.new)
  
  v2.new <- B1.p %*% v2.old + B2.p %*% v2.old
  v2.new <- v2.new - as.numeric( (t(v1.new) %*% v2.new) / (t(v1.new) %*% v1.new) ) * v1.new
  v2.new <- v2.new/norm(v2.new)
  
  i <- 1
  
  while( i < 100) {
    v1.old <- v1.new
    v2.old <- v2.new
    
    v1.new <- B1.p %*% v1.old + B2.p %*% v1.old
    v1.new <- v1.new / norm(v1.new)
    
    v2.new <- B1.p %*% v2.old + B2.p %*% v2.old
    v2.new <- v2.new <- v2.new - as.numeric( (t(v1.new) %*% v2.new) / (t(v1.new) %*% v1.new) ) * v1.new
    v2.new <- v2.new/norm(v2.new)
    i <- i + 1
  }
  return(list(vector1 = v1.new, 
              value1 = t(v1.new) %*% (B1.p + B2.p) %*% v1.new, 
              topnode1 = which.max(v1.new),
              vector2 = v2.new, 
              value2 = t(v2.new) %*% (B1.p + B2.p) %*% v2.new, 
              topnode2 = which.max(v2.new),
              iter = i
  ))
}

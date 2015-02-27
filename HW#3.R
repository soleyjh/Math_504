# JAMES HYER SOLEY
# NUMERICAL METHODS
# MATH 504
# 02/01/2015
# HOMEWORK #3

# Problem #2
options(digits=22)
x <- 1/3
x
x <- x + 1
x <- x - 1
x
x == 1/3
print(x)

# Part A
x <- 1/3
x
x <- x + 10^4
x
x <- x - 10^4
x
x == 1/3
print(x)

# Part B
x <- 1/3
x
x <- x + 10^15
x
x <- x - 10^15
x
x == 1/3
print(x)

# Problem #4
# Part A
proj.u <- function(u, v){
  c <- solve(t(u)%*%U, t(u)%*%v)
  u.star <- U%*%c
  return(u.star)
}

# Part B
e <- function(v,u){norm(v-proj(u,v))}

# Problem #5
lines.dist <- function(x.1,x.2,y.1,y.2){
  s1 <- seq(-1000,1000,1)
  s2 <- seq(-1000,1000,1)
  
  # In 2-d the functions must be coplanar and thus cannot be skew lines
  # therefore they either intersect or are parrallel by definition
  
  if (length(x.1) < 3){
    if (identical(x.2,y.2)){
      # Parrallel Case
      s.1 <- sample(s1,1)
      s.2 <- sample(s2,1)
      d <- sqrt((x.2 - x.1)^2 + (y.2 - y.1)^2)
      return(d)
      print("Lines are Parallel in 2D so s1 and s2 can be anything")
    }
    else{
      # Intersect Case 
      d <- 0
      return(d)
      print("Lines are NOT Parallel in 2D so lines must Intersect and D = 0")
    }
  }
  # In GT 2-d the functions are skew lines
  # Therefore for range of S values use distance formula in R = d
    
  else{
    dist <- function(s1, s2) {   
    d <- sqrt(t(x.1 + s1*x.2 - (y.1 + s2*y.2)) %*% (x.1 + s1*x.2 - (y.1 + s2*y.2))) 
    return(d)}
    
    start <- cbind(seq(-1000,1000,1),seq(-1000,1000,1)) 
    m<-min(mapply(function(s1,s2) nlm(dist,s1,s2)$minimum,start[,1],start[,2]))
    return(m)}  
}

#3D
x.1 <- c(-2,1,0)
x.2 <- c(1,-1,1)
y.1 <- c(0,1,0)
y.2 <- c(1,1,2)

lines.dist(x.1,x.2,y.1,y.2)

#2D Intersect
x.1 <- c(-2,1)
x.2 <- c(1,-1)
y.1 <- c(0,1)
y.2 <- c(1,1)

lines.dist(x.1,x.2,y.1,y.2)

#2d Parallel
x.1 <- c(-2,1)
x.2 <- c(1,-1)
y.1 <- c(0,1)
y.2 <- c(1,-1)
lines.dist(x.1,x.2,y.1,y.2)
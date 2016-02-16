# JAMES HYER SOLEY
# NUMERICAL METHODS
# MATH 611
# 03/18/2015

# Problem #2
# Part A

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



le <- read.delim("~/Math_504/HW#7/life_expectancy.txt")
x <- c(41:100);
plot(le[x,2],cex=.5,col="darkgrey", main='Spline Fit', ylab="Probability a person Dies", xlab="Age")
l <- spline(le[x,2])
lines(l)

# Part B

y <- cumprod(1-le[x,2])
i <- c(1:720);
L <- splinefun(x,y)(40+i/12)
sum(200*L*exp(-0.05*i/12))

# Problem #3
# Part BiiA

x1 <- 15
x2 <- 20

r1 <- c(1,x1,x1^2,x1^3,-1,-x1,-x1^2,-x1^3,0,0,0,0)
r2 <- c(0,1,2*x1,3*x1^2,0,-1,-2*x1,-3*x1^2,0,0,0,0)
r3 <- c(0,0,2,6*x1,0,0,-2,-6*x1,0,0,0,0)
r4 <- c(0,0,0,0,1,x2,x2^2,x2^3,-1,-x2,-x2^2,-x2^3)
r5 <- c(0,0,0,0,0,1,2*x2,3*x2^2,0,-1,-2*x2,-3*x2^2)
r6 <- c(0,0,0,0,0,0,2,6*x2,0,0,-2,-6*x2)

M <- rbind(r1,r2,r3,r4,r5,r6)
M

z.bar <- rep(1,12)
M %*% z.bar # equals 0 vector, as expected
M.t <- t(M)
M.t
z.bar2 <- rep(1,6)
M.t %*% z.bar2

# Part BiiB
bones <- read.table("~/Math_504/HW#7/bones.txt", header=TRUE, quote="\"")
male <- bones[which(bones$gender == 'male'), ]
female <- bones[which(bones$gender == 'female'), ]

M.t <- t(M)
A <- cbind(M.t, matrix(rnorm(72), nrow=12, ncol=6))
det(A)

A.qr <- qr(A)
Q <- qr.Q(A.qr)
Q <- Q[,7:12]

g.x <- function(z = Q[,1], data = male$age) {
  store <- NULL
  for(i in 1:length(data)){
    if(data[i] < 15) {
      s.x1 <- z[1] + z[2]*data[i] + z[3]*data[i]^2 + z[4]*data[i]^3
      store <- c(store, s.x1)
    } else if (data[i] < 20 & data[i] >=15) {
      s.x2 <- z[5] + z[6]*data[i] + z[7]*data[i]^2 + z[8]*data[i]^3
      store <- c(store, s.x2)
    } else if (data[i] >= 20) {
      s.x3 <- z[9] + z[10]*data[i] + z[11]*data[i]^2 + z[12]*data[i]^3
      store <- c(store, s.x3)
    }
  }
  return(store)
}

fit.male <- lapply(1:6, function(x) g.x(z = Q[,x], data = male$age))

# Part C
library(ggplot2)
h <- list(
  h.1 = function(x) return(1),
  h.2 = function(x) return(x),
  h.3 = function(x) return(x^2),
  h.4 = function(x) return(x^3),
  h.5 = function(x) return(ifelse((x - 15)^3 > 0, (x - 15)^3, 0)),
  h.6 = function(x) return(ifelse((x - 20)^3 > 0, (x - 20)^3, 0))
)

h.i <- lapply(1:6, function(y) h[[y]](x=male$age))
h.plots <- lapply(1:6, function(x) qplot(male$age,unlist(h.i[[x]]), geom="line") + ylab(x))
g.plots <- lapply(1:6, function(x) qplot(male$age,unlist(fit.male[[x]]), geom="line") + ylab(x))
plots <- append(h.plots, g.plots)
multiplot(plotlist = plots, layout=matrix(c(1:12), nrow=6, ncol=2))

# Part D
y.male <- male$spnbmd
h.i <- lapply(1:6, function(y) h[[y]](x=male$age))
W.male <- cbind(unlist(h.i[[1]]), unlist(h.i[[2]]), unlist(h.i[[3]]), unlist(h.i[[4]]), unlist(h.i[[5]]), unlist(h.i[[6]])) 

alpha.male <- solve(t(W.male)%*%W.male) %*% t(W.male) %*% y.male 

fit.hm <- W.male %*% alpha.male
qplot(male$age, y.male) + geom_line(aes(y=fit.hm, color="red"))+ theme(legend.position="none")

y.female <- female$spnbmd
h.i <- lapply(1:6, function(y) h[[y]](x=female$age))
W.female <- cbind(unlist(h.i[[1]]), unlist(h.i[[2]]), unlist(h.i[[3]]), unlist(h.i[[4]]), unlist(h.i[[5]]), unlist(h.i[[6]])) 

alpha.female <- solve(t(W.female)%*%W.female) %*% t(W.female) %*% y.female 

fit.hf <- W.female %*% alpha.female
qplot(female$age, y.female) + geom_line(aes(y=fit.hf, color=3))+ theme(legend.position="none")

# Part E
y.male <- male$spnbmd
g.i <- lapply(1:6, function(x) g.x(z = Q[,x], data = male$age))
W.gm <- cbind(unlist(g.i[[1]]), unlist(g.i[[2]]), unlist(g.i[[3]]), unlist(g.i[[4]]), unlist(g.i[[5]]), unlist(g.i[[6]])) 
alpha.gm <- solve(t(W.gm)%*%W.gm) %*% t(W.gm) %*% y.male 
fit.gm <- W.gm %*% alpha.gm

qplot(male$age, y.male) + geom_line(aes(y=fit.gm, color="red", size=.3)) + geom_line(aes(y=fit.hm, color="green", size = .1)) + theme(legend.position="none")


y.female <- female$spnbmd
g.i <- lapply(1:6, function(x) g.x(z = Q[,x], data = female$age))
W.gf <- cbind(unlist(g.i[[1]]), unlist(g.i[[2]]), unlist(g.i[[3]]), unlist(g.i[[4]]), unlist(g.i[[5]]), unlist(g.i[[6]])) 
alpha.gf <- solve(t(W.gf) %*% W.gf) %*% t(W.gf) %*% y.female 

fit.gf <- W.gf %*% alpha.gf

qplot(female$age, y.female) + geom_line(aes(y=fit.gf, color="blue", size=.1)) + geom_line(aes(y=fit.hf, size=.01)) + theme(legend.position="none")

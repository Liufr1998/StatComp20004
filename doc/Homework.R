## -----------------------------------------------------------------------------
boxplot(count ~ spray, data = InsectSprays, col = "lightgray")

## -----------------------------------------------------------------------------
library(DT)
datatable(mtcars)

## -----------------------------------------------------------------------------
set.seed(100000)
a <- 2
b <- 2
n <- 1000
u <- runif(n)
x <- b/((1-u)^(1/a))
hist(x,prob=T,breaks = seq(2,ceiling(max(x)/10)*10+2,by=10),main=expression(f(x)==frac(ab^a,x^(a+1))))
y <- seq(2,max(x),length=100)
lines(y,a*(b^a)/(y^(a+1)))

## -----------------------------------------------------------------------------
f1 <- function(n){
  u1 <- runif(n, min=-1, max=1)
  u2 <- runif(n, min=-1, max=1)
  u3 <- runif(n, min=-1, max=1)
  u <- numeric()
  for(i in 1:n)
    if(abs(u3[i])>=max(abs(u2[i]),abs(u1[i]))) u[i] <- u2[i] else u[i] <- u3[i]
  return(u)
}
n <- 100000
x <- f1(n)
hist(x,pro=T,ylim=c(0,0.8),main=expression(f[e](x)==frac(3,4)(1-x^2)))
lines(density(x,kernel="epanechnikov"),col="blue",lwd=3)
y <- seq(-1,1,by=0.05)
lines(y,3/4*(1-y^2),col="red",lwd=3)

## -----------------------------------------------------------------------------
set.seed(100000)
n <- 1000
r <- 4
beta <- 2
lambda <- rgamma(n,r,beta)
x <- rexp(n,lambda)
hist(x,pro=T,ylim=c(0,1.4),main=expression(f(x)==frac(r*beta^r,(beta+x)^(r+1))))
lines(density(x,kernel="gaussian"),col="blue",lwd=3)
y <- seq(0,max(x),length=100)
lines(y,r*beta^r/(beta+y)^(r+1),col="red",lwd=3)

## -----------------------------------------------------------------------------
set.seed(12345)
ns = 10000
u = runif(ns,0,pi/3)
estimator = (pi/3)*mean(sin(u))
exact = 0.5
print(c(estimator, exact, estimator-exact))

## -----------------------------------------------------------------------------
set.seed(12345)
ns = 10000
nm = 1000
theta1 = theta2 = numeric(nm)
for(ii in 1:nm){
  u1 = runif(ns,0,1)
  u2 = runif(ns/2,0,1)
  theta1[ii] = mean(exp(u1)) # simple Monte Carlo
  theta2[ii] = mean(exp(u2)+exp(1-u2))/2 # antithetic variate approach
}
var1 = var(theta1)
var2 = var(theta2)
p_estimated = (var1-var2)/var1
p_thero = 2*(exp(2)-3*exp(1)+1)/(-exp(2)+4*exp(1)-3)
print(c(p_estimated,p_thero))

## -----------------------------------------------------------------------------
set.seed(12345)
M <- 1e4
g <- function(x) x^2*exp(-x^2/2)/sqrt(2*pi)*(x>1)
theta.hat <- se <- numeric(2)

u <- runif(M)
x1 <- sqrt(-2*log(1-u)+1)
fg1 <- g(x1)/(x1*exp(-(x1^2-1)/2))*(x1>1)
theta.hat[1]<- mean(fg1)
se[1] <- sd(fg1)

x2 <- 1-log(1-u)
fg2 <- g(x2)/exp(-x2+1)*(x2>1)
theta.hat[2]<- mean(fg2)
se[2] <- sd(fg2)

matrix(c(theta.hat,se), nrow = 2, dimnames = list(c("f1","f2"),c("theta.hat","se")))

## -----------------------------------------------------------------------------
set.seed(12345)
M <- 10000
k <- 5 # the number of stratum
r <- M/k # repetitions in each stratum
N <- 50 # the number of times to repeat the estimation
estimator <- matrix(nrow=N, ncol=k)
f1 <- function(x,j) (exp(-(j-1)/k)-exp(-j/k))/(1+x^2)*(x>(j-1)/k)*(x<j/k)
for(i in 1:N)
  for(j in 1:k)
    estimator[i,j] <- mean(f1(runif(r,(j-1)/k,j/k),j))
est <- sum(apply(estimator,2,mean))
est.var <- sum(apply(estimator,2,var))/r
est.sd <- sqrt(est.var)
est.sd1 <- 0.09658794/sqrt(M)
print(c(0.52506988, est))
print(c(est.sd1, est.sd))

## -----------------------------------------------------------------------------
set.seed(12345)
mu <- 0; sig2 <- 1 # parameters
m <- 1e4 # repetitions
n <- 50 # sample size
alpha <- 0.05
mu.hat <- matrix(nrow=m, ncol=2)
T1 <- numeric(m)
for(i in 1:m){
  lx <- log(rlnorm(n, mu, sig2))
  mu.hat[i,1] <- mean(lx)+qt(alpha/2,n-1)*sd(lx)/sqrt(n)
  mu.hat[i,2] <- mean(lx)-qt(alpha/2,n-1)*sd(lx)/sqrt(n)
  T1[i] <- (mu.hat[i,1]<=mu)*(mu.hat[i,2]>=mu)
}
mean(T1)

## -----------------------------------------------------------------------------
set.seed(12345)
mu <- 2; sig2 <- 2*mu # parameters
m <- 1e4 # repetitions
n <- 20 # sample size
alpha <- 0.05
mu.hat <- matrix(nrow=m, ncol=2)
sig2.hat <- T1 <- T2 <- numeric()
for(i in 1:m){
  x <- rchisq(n,mu)
  mu.hat[i,1] <- mean(x)+qt(alpha/2,n-1)*sd(x)/sqrt(n)
  mu.hat[i,2] <- mean(x)-qt(alpha/2,n-1)*sd(x)/sqrt(n)
  T1[i] <- (mu.hat[i,1]<=mu)*(mu.hat[i,2]>=mu)
  sig2.hat <- (n-1)*var(x)/qchisq(alpha, n-1)
  T2[i] <- sig2.hat>=sig2
}
mean(T1)
mean(T2)

## ---- eval=FALSE--------------------------------------------------------------
#  set.seed(12345)
#  alpha <- 0.05 # significant level
#  n <- 100 # sample size
#  m <- 1000 # repetition times
#  cv <- qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3)))) # critical value
#  pow.b <- pow.t <- numeric()
#  
#  sk <- function(x){
#    xbar <- mean(x)
#    mean((x-xbar)^3)/(mean((x-xbar)^2))^1.5
#  }
#  
#  theta <- seq(1,50,by=1) # parameter
#  for(i in 1:length(theta)){ # for each paramater
#    sktest.b <- sktest.t <- numeric()
#    for(j in 1:m){ # for each replicate
#      x <- rbeta(n,theta[i],theta[i])
#      sktest.b[j] <- as.integer(abs(sk(x))>=cv)
#      x <- rt(n,theta[i])
#      sktest.t[j] <- as.integer(abs(sk(x))>=cv)
#    }
#    pow.b[i] <- mean(sktest.b)
#    pow.t[i] <- mean(sktest.t)
#  }
#  
#  plot(theta,pow.b,type="b",xlab=bquote(alpha),ylab="power")
#  abline(h=alpha,lty=3)
#  plot(theta,pow.t,type="b",xlab=bquote(nu),ylab="power")
#  abline(h=alpha,lty=3)

## ---- eval=FALSE--------------------------------------------------------------
#  count5test <- function(x, y) {
#    X <- x - mean(x)
#    Y <- y - mean(y)
#    outx <- sum(X > max(Y)) + sum(X < min(Y))
#    outy <- sum(Y > max(X)) + sum(Y < min(X))
#    return(as.integer(max(c(outx, outy)) > 5))
#  }
#  
#  Ftest <- function(x,y) {
#    alpha <- 0.055
#    f <- var(x)/var(y)
#    m1 <- length(x)
#    n1 <- length(y)
#    return(as.integer(f<=pf(alpha/2,m1-1,n1-1) | f>=pf(1-alpha/2,m1-1,n1-1)))
#  }
#  
#  set.seed(12345)
#  m <- 10000
#  n <- c(30,100,500)
#  sigma1 <- 1
#  sigma2 <- 1.5
#  power <- matrix(nrow=2,ncol=3,dimnames=list(c("Count Five test","F test"),c("n=30","n=100","n=500")))
#  
#  for(i in 1:length(n)){
#    out.c <- out.f <- numeric()
#    for(j in 1:m){
#      x <- rnorm(n[i], 0, sigma1)
#      y <- rnorm(n[i], 0, sigma2)
#      out.c[j] <- count5test(x, y)
#      out.f[j] <- Ftest(x,y)
#    }
#  
#    power[1,i] <- mean(out.c)
#    power[2,i] <- mean(out.f)
#  }
#  print(power)

## ---- eval=FALSE--------------------------------------------------------------
#  multi.sk <- function(x){
#    n <- nrow(x)
#    xbar <- colMeans(x)
#    sig.inv <- solve(cov(x)*(n-1)/n)
#    b <- 0
#    for (i in 1:n) {
#      for (j in 1:n) {
#        b <- b+(t(x[i,]-xbar)%*%sig.inv%*%(x[j,]-xbar))^3
#      }
#    }
#    return(b/(n^2))
#  }
#  
#  #example 6.8
#  library(MASS)
#  set.seed(12345)
#  n <- c(10,20,30,50,100,500)
#  m <- 1000
#  d <- 2
#  sigma <- matrix(c(1,0,0,1),nrow=d)
#  alpha <- 0.05
#  cv <- cbind((6/n)*qchisq(alpha/2,d*(d+1)*(d+2)/6),
#              (6/n)*qchisq((1-alpha/2),d*(d+1)*(d+2)/6))
#  p <- numeric()
#  for (i in 1:length(n)) {
#    msktests <- numeric()
#    for (j in 1:m) {
#      x <- mvrnorm(n[i],rep(0,2),sigma)
#      t <- multi.sk(x)
#      msktests[j] <- as.integer(t<=cv[i,1] | t>=cv[i,2])
#    }
#    p[i]<- mean(msktests)
#  }
#  matrix(c(n,p),nrow=2,dimnames=list(c("n","estimate"),1:length(n)),byrow=T)
#  
#  #example 6.10
#  n <- 30
#  m <- 2500
#  alpha <- 0.1
#  sigma <- list(matrix(c(1,0,0,1),d),
#                matrix(c(10,0,0,10),d))
#  epsilon <- c(seq(0,0.15,0.01),seq(0.15,1,0.05))
#  N <- length(epsilon)
#  mpwr <- numeric()
#  cv <- c(alpha/2,d*(d+1)*(d+2)/6,(6/n)*qchisq((1-alpha/2),d*(d+1)*(d+2)/6))
#  for (i in 1:N) {
#    e <- epsilon[i]
#    msktests <- numeric()
#    for (j in 1:m) {
#      k <- sample(c(1,2),replace = TRUE,size = n,prob = c(e,1-e))
#      y <- matrix(nrow=n,ncol=d)
#      x <- list(mvrnorm(n, rep(0,2), sigma[[1]]),
#                mvrnorm(n, rep(0,2), sigma[[2]]))
#      for (p in 1:n) y[p,] <- x[[k[p]]][p,]
#      t <- multi.sk(y)
#      msktests[j] <- as.integer(t<=cv[1] | t>=cv[2])
#    }
#    mpwr[i]<- mean(msktests)
#  }
#  #plot power vs epsilon
#  plot(epsilon, mpwr, type = "b",
#       xlab = bquote(epsilon), ylim = c(0,1))
#  abline(h = .1, lty = 3)
#  se <- sqrt(mpwr * (1-mpwr) / m) #add standard errors
#  lines(epsilon, mpwr+se, lty = 3)
#  lines(epsilon, mpwr-se, lty = 3)

## ---- eval=FALSE--------------------------------------------------------------
#  S<- (0.651*(1-0.651)+(1-0.676)*0.676)/10000
#  c(qnorm(0.025,0.025,S),qnorm(0.975,0.025,S))

## -----------------------------------------------------------------------------
library("bootstrap")
data(law)
n <- nrow(law)
cor.hat <- cor(law$LSAT,law$GPA)
est <- numeric()
for(i in 1:n) est[i] <- cor(law$LSAT[-i],law$GPA[-i])
bias <- (n-1)*(mean(est)-cor.hat)
se <- sqrt((n-1)*mean((est-mean(est))^2))
matrix(c(bias,se),nrow=1,dimnames=list("Jackknife estimate",c("bias","standard error")))

## -----------------------------------------------------------------------------
set.seed(12345)
library("boot")
data("aircondit")
boot.mean <- function(x,i) mean(x[i])
ci.norm <- ci.basic <- ci.perc <- ci.bca <- numeric()
de <- boot(data=unlist(aircondit), statistic=boot.mean, R=2000)
ci <- boot.ci(de, type=c("norm","basic","perc","bca"))
matrix(c(ci$norm[2:3], ci$basic[4:5], ci$percent[4:5], ci$bca[4:5]), ncol=2, byrow=T,
       dimnames=list(c("standard normal", "basic", "percentile", "BCa"), c("Lower limit", "Upper limit")))


## -----------------------------------------------------------------------------
library("bootstrap")
data("scor")
theta.hat <- eigen(cov(scor))$value[1]/sum(diag(cov(scor)))
n <- nrow(scor)
theta.J <- numeric(n)
for(i in 1:n){
  x <- scor[-i,]
  theta.J[i] <- eigen(cov(x))$value[1]/sum(diag(cov(x))) 
}
bias.J <- (n-1)*(mean(theta.J)-theta.hat)
se.J <- sqrt((n-1)*mean((theta.J-mean(theta.J))^2))
matrix(c(bias.J,se.J), nrow = 1,
       dimnames=list("Jackknife estimate",c("bias","standard error")))

## -----------------------------------------------------------------------------
library(DAAG); attach(ironslag)
n <- length(magnetic)
e <- rep(0,4)
# for n-fold cross validation
# fit models on leave-two-out samples
for (i in 1:(n-1)) {
  for(j in (i+1):n){
    
    y <- magnetic[-c(i,j)]
    x <- chemical[-c(i,j)]
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[c(i,j)]
    e[1] <- e[1] + sum((magnetic[c(i,j)] - yhat1)^2)
    
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[c(i,j)] +
      J2$coef[3] * chemical[c(i,j)]^2
    e[2] <- e[2] + sum((magnetic[c(i,j)] - yhat2)^2)
    
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(i,j)]
    yhat3 <- exp(logyhat3)
    e[3] <- e[3] + sum((magnetic[c(i,j)] - yhat3)^2)
    
    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(i,j)])
    yhat4 <- exp(logyhat4)
    e[4] <- e[4] + sum((magnetic[c(i,j)] - yhat4)^2)
    
  }
}
matrix(e/n/(n-1), nrow=1,
       dimnames=list("prediction error", c("Linear","Quadratic"," Exponential","Log-Log")))


## ---- eval=FALSE--------------------------------------------------------------
#  count5test <- function(x, y) {
#    X <- x - mean(x)
#    Y <- y - mean(y)
#    outx <- sum(X > max(Y)) + sum(X < min(Y))
#    outy <- sum(Y > max(X)) + sum(Y < min(X))
#    return(max(c(outx, outy)))
#  }
#  
#  test <- function(x,y){
#    z <- c(x,y)
#    m <- length(x)
#    N <- length(z)
#    D <- numeric()
#    D0 <- count5test(x,y)
#  
#    for(i in 1:R){
#      k <- sample(1:N, m, replace = F)
#      x1 <- z[k]
#      y1 <- z[-k]
#      D[i] <- count5test(x1,y1)
#    }
#  
#    as.integer(mean(c(D0,D)>D0) < alpha)
#  }
#  
#  
#  n1 <- 20
#  n2 <- 30
#  mu1 <- mu2 <- 0
#  sigma1 <- sigma2 <- 1
#  m <- 10000
#  R <- 999
#  alpha <- 0.05
#  set.seed(12345)
#  alphahat <- mean(replicate(m,expr={
#    x <- rnorm(n1, mu1, sigma1)
#    y <- rnorm(n2, mu2, sigma2)
#    test(x,y)
#  }))
#  alphahat
#  

## ---- eval=FALSE--------------------------------------------------------------
#  library(boot)
#  library(RANN)
#  library(energy)
#  library(Ball)
#  
#  Tn <- function(z, ix, sizes,k) {
#    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z,0);
#    z <- z[ix, ];
#    NN <- nn2(data=z, k=k+1) # what's the first column?
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#    (i1 + i2) / (k * n)
#  }
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=R,
#                     sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)
#  }

## ---- eval=FALSE--------------------------------------------------------------
#  # Unequal variances and equal expectations
#  m <- 1e3; k<-2; p<-2; sigma1<-1; sigma2<-1.5; set.seed(12345)
#  n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p,sd=sigma1),ncol=p)
#    y <- matrix(rnorm(n2*p,sd=sigma2),ncol=p)
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=R,seed=i*12345)$p.value
#  }
#  alpha <- 0.1
#  pow <- colMeans(p.values<alpha)
#  matrix(pow,nrow=1,dimnames=list("power",c("NN","energy","Ball")))

## ---- eval=FALSE--------------------------------------------------------------
#  # Unequal variances and unequal expectations
#  m <- 1e3; k<-1; p<-1; mu1<-0; mu2<-0.5; sigma1<-1; sigma2<-1.5; set.seed(12345)
#  n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p,mu1,sigma1),ncol=p)
#    y <- matrix(rnorm(n2*p,mu2,sigma2),ncol=p)
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=R,seed=i*12345)$p.value
#  }
#  alpha <- 0.05
#  pow <- colMeans(p.values<alpha)
#  matrix(pow,nrow=1,dimnames=list("power",c("NN","energy","Ball")))

## ---- eval=FALSE--------------------------------------------------------------
#  # Non-normal distributions
#  m <- 1e3; k<-1; p<-1; mu1<-0; mu2<-0.5; sigma1<-1; sigma2<-1; epi<-0.5; set.seed(12345)
#  n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rt(n1*p,1),ncol=p)
#    y <- matrix(epi*rnorm(n2*p,mu1,sigma1)+(1-epi)*rnorm(n2*p,mu2,sigma2),ncol=p)
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=R,seed=i*12345)$p.value
#  }
#  alpha <- 0.05
#  pow <- colMeans(p.values<alpha)
#  matrix(pow,nrow=1,dimnames=list("power",c("NN","energy","Ball")))

## ---- eval=FALSE--------------------------------------------------------------
#  # Unbalanced samples
#  m <- 1e3; k<-2; p<-2; sigma1<-1; sigma2<-1.5; set.seed(12345)
#  n1<-50; n2<-500; R<-999; n <- n1+n2; N = c(n1,n2)
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p,sd=sigma1),ncol=p)
#    y <- matrix(rnorm(n2*p,sd=sigma2),ncol=p)
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=R,seed=i*12345)$p.value
#  }
#  alpha <- 0.1
#  pow <- colMeans(p.values<alpha)
#  matrix(pow,nrow=1,dimnames=list("power",c("NN","energy","Ball")))
#  

## -----------------------------------------------------------------------------
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (exp(-abs(y))/2) / (exp(-abs(x[i-1]))/2))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      } }
  return(list(x=x, k=k))
}

set.seed(12345)
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)

print(c(rw1$k, rw2$k, rw3$k, rw4$k)/N)

a <- c(.025, .975)
rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
mc <- rw[501:N, ]
Qrw <- apply(mc, 2, function(x) quantile(x, a))

# par(mfrow=c(2,2))  #display 4 graphs together
rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
for (j in 1:4) {
  plot(rw[,j], type="l",
       xlab=bquote(sigma == .(round(sigma[j],3))),
       ylab="X", ylim=range(rw[,j]))
  abline(h=Qrw[,j])
}
par(mfrow=c(1,1)) #reset to default



## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  
  psi.means <- rowMeans(psi)     #row means
  B <- n * var(psi.means)        #between variance est.
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within est.
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
  r.hat <- v.hat / W             #G-R statistic
  return(r.hat)
}

Laplace.chain <- function(sigma, N, X1) {
  x <- rep(0, N)
  x[1] <- X1
  u <- runif(N)
  
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (exp(-abs(y))/2) / (exp(-abs(x[i-1]))/2))
      x[i] <- y else x[i] <- x[i-1]
  }
  return(x)
}

sigma <- 1     #parameter of proposal distribution
k <- 4          #number of chains to generate
n <- 15000      #length of chains
b <- 1000       #burn-in length

#choose overdispersed initial values
x0 <- c(-10, -5, 5, 10)

#generate the chains
set.seed(12345)
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
  X[i, ] <- Laplace.chain(sigma, n, x0[i])

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))

#plot psi for the four chains
#    par(mfrow=c(2,2))
for (i in 1:k)
  if(i==1){
    plot((b+1):n,psi[i, (b+1):n],ylim=c(-0.2,0.2), type="l",
         xlab='Index', ylab=bquote(phi))
  }else{
    lines(psi[i, (b+1):n], col=i)
  }
par(mfrow=c(1,1)) #restore default


#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
kopt <- c(4:25,100,500,1000)
Ak <- numeric()
set.seed(12345)
for(i in 1:length(kopt)){
  k <- kopt[i]
  f <- function(a) pt(sqrt(a^2*(k-1)/(k-a^2)),df=k-1)-pt(sqrt(a^2*k/(k+1-a^2)),df=k)
  res <- uniroot(f,c(0,sqrt(k))-1e-10,maxiter = 1e6)
  res <- as.numeric(unlist(res)[1])
  if(res>0 & res<sqrt(k)) Ak[i] <- res else Ak[i] <- NA
}
matrix(Ak,ncol=1,dimnames = list(paste("k=",kopt,sep=""),"A(k)"))

## -----------------------------------------------------------------------------
EL <- function(pq)
  -(((po^2/(po^2+2*po*(1-po-qo))+1)*nA.+nAB)*log(pq[1])+
   ((qo^2/(qo^2+2*qo*(1-po-qo))+1)*nB.+nAB)*log(pq[2])+
   ((1-po^2/(po^2+2*po*(1-po-qo)))*nA.+(1-qo^2/(qo^2+2*qo*(1-po-qo)))*nB.+2*nOO)*log(1-sum(pq)))
   
nA.<-444;nB.<-132;nOO<-361;nAB<-63
iter <- 20
po <- qo <- 0.2
iter.result <- matrix(nrow=iter+1,ncol=2,
                      dimnames=list(c("initial",paste("iteration",1:iter,sep="")),c("p","q")))
logml <- numeric()
iter.result[1,] <- c(po,qo)
logml[1] <- -EL(c(po,qo))

for(i in 1:iter){
  if(i>1) {po <- iter.result[i-1,1]; qo <- iter.result[i-1,2]}
  em <- optim(c(po,qo), EL)
  iter.result[i+1,] <- em$par
  logml[i+1] <- -em$value
}
cbind(iter.result,logml)


## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

# loops
n <- length(formulas)
result <- list()
for(i in 1:n) result[[i]] <- lm(formulas[[i]], data = mtcars)
result

# lapply
lapply(formulas, lm, data = mtcars)

## -----------------------------------------------------------------------------
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
) 
sapply(trials,function(x) x$p.value)

## -----------------------------------------------------------------------------
f <- function(data,fun,FUN.VALUE){
  data2 <- Map(fun,data)
  vapply(data2,function(x) x,FUN.VALUE)
}

FUN.VALUE <- double(1)
fun <- mean
data <- data.frame(replicate(6,sample(c(1:10),10,rep=T)))
f(data,fun,FUN.VALUE)

## -----------------------------------------------------------------------------
library(Rcpp)
# sourceCpp("CMetro.cpp") 

set.seed(3000)

lap_f = function(x) exp(-abs(x))

RrwM = function(sigma, x0, N){
  x = numeric(N)
  x[1] = x0
  u = runif(N)
  k = 0
  for (i in 2:N) {
    y = rnorm(1, x[i-1], sigma)
    if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) x[i] = y 
    else {
      x[i] = x[i-1]
      k = k+1
    }
  }
  return(list(x = x, k = k))
}

N = 2000
sigma = c(.05, .5, 2, 16)
x0 = 25
crw = list()
rrw = list()
for(i in 1:length(sigma)){
  crw[[i]] = CrwM(sigma[i],x0,N)
  rrw[[i]] = RrwM(sigma[i],x0,N)
}

# R function
# par(mfrow=c(2,2)) 
for (j in 1:4) {
  plot(rrw[[j]]$x, type="l",
       xlab=bquote(sigma == .(round(sigma[j],3))),
       ylab="X", ylim=range(rrw[[j]]$x))
}

# Rcpp function
# par(mfrow=c(2,2))
for (j in 1:4) {
  plot(crw[[j]], type="l",
       xlab=bquote(sigma == .(round(sigma[j],3))),
       ylab="X", ylim=range(crw[[j]]))
}

# Q-Q plot
# par(mfrow=c(2,2))
for(i in 1:length(sigma)){
  qqplot(rrw[[i]]$x,crw[[i]],xlab="RrwM",ylab="CrwM",
         main=bquote(sigma == .(round(sigma[i],3))))
  abline(0,1,col="red")
}

## -----------------------------------------------------------------------------
library(microbenchmark)
microbenchmark(rwR=RrwM(sigma[1],x0,N),rwC=CrwM(sigma[1],x0,N))


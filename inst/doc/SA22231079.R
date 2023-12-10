## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
mysample <- function(str=1,stp=100, num = 20){ #默认从1到100中有放回地随机挑选20个变量
  n = stp - str + 1
  x = seq(str, stp, 1)
  p = rep(1/n, n)
  Fu = c(0,cumsum(p));  #累加得到分布函数
  y = integer(num)  
  m = length(Fu)
  u = runif(num)  
  for(i in 1:m-1){   #逆转函数
    index = u>Fu[i] & u<=Fu[i+1]; 
    y[index]=x[i];
  };
  y
}

## -----------------------------------------------------------------------------
U = runif(1000,0,1)
x = integer(1000)
for(i in 1:1000){      #reverse func
  if(U[i] > 0.5) {
    x[i] = -log(2-2*U[i]) 
    }
  else if(U[i] < 0.5) {
    x[i] = log(2*U[i])
      }
  else U[i] = 0.5
}
hist(x, prob = TRUE,ylim = c(0,0.8))  # Plot the frequency histogram for generating a random number x
xx = seq(-100,100,0.1)
lines(xx ,0.5 * exp(-abs(xx))) #Curve Fitting Controls

## -----------------------------------------------------------------------------
myB <- function(n=1000,a=3,b=2){      #Default generation of 1000 random variables conforming to B(3, 2) distribution
  F <- function(x) (x^(a-1))*((1-x)^(b-1))
  Be = integrate(F,0,1)
  i=0 
  x = integer(n)
  while(i<n)
  {
    u = runif(1)
    y = runif(1)
    if(F(y) > u){                     #A-C
      i = i+1
      x[i] = y
    }
  }
  x
}
B <- myB()
hist(B,prob=T)  
xx = seq(0,1,0.01)
lines(xx,dbeta(xx,3,2)) # Curve Fitting Controls

## -----------------------------------------------------------------------------
F <- function(x) 0.75*(1-x^2)   #topic-specific recapitulation
n <- 1000  
u <- integer(n)
for (i in 1:n) {
  u1 <- runif(n,-1,1)
  u2 <- runif(n,-1,1)
  u3 <- runif(n,-1,1)
  ifelse(abs(u3[i]) >= abs(u2[i]) && abs(u3[i]) >= abs(u1[i]),u[i] <- u2[i],
         u[i] <- u3[i])
}
hist(u, prob = TRUE)
x <- seq(-1, 1, 0.01)
lines(x, F(x))

## -----------------------------------------------------------------------------
n <- 1E6  
K <- 100  

rho_1 <- 1.0  # ρmin
rho_2 <- 0.5  # ρ2
rho_3 <- 0.8  # ρ3

# simulate pi
simulate_pi <- function(rho) {
  pi_values <- numeric(K)
  for (i in 1:K) {
    m <- rbinom(1, n, rho)
    # calcu p_hat and pi_hat
    p_hat <- n / m
    pi_hat <- 2 * rho * p_hat
    pi_values[i] <- pi_hat
  }
  return(pi_values)
}


pi_values_1 <- simulate_pi(rho_1)
var_pi_1 <- var(pi_values_1)

pi_values_2 <- simulate_pi(rho_2)
var_pi_2 <- var(pi_values_2)

pi_values_3 <- simulate_pi(rho_3)
var_pi_3 <- var(pi_values_3)


print(paste("ρmin：", var_pi_1))
print(paste("ρ2:", var_pi_2))
print(paste("ρ3:", var_pi_3))

## -----------------------------------------------------------------------------
m <- 1E4
simula <- function(){
  x <- runif(m, min=0, max=1)
  theta1.hat <<- mean(exp(x)) 
  x1 <- runif(m/2, min=0, max=1)
  theta2.hat <<- (mean(exp(x1))+mean(exp(1-x1)))/2 
  true <- exp(1)-1
  c(theta1.hat,theta2.hat,true)
}
simula()


V1 = integer(m)
V2 = integer(m)
for (i in 1:m){
  simula()
  V1[i] <- theta1.hat
  V2[i] <- theta2.hat
}

a <- var(V1)
b <- var(V2)
ans <- (a-b)/a

## -----------------------------------------------------------------------------
m <- 1e4
f1 <- function(u)u^(1/4)
f2 <- function(u)exp(-(1/2*u^2))
set.seed(100)
u <- runif(m)
T1 <- f1(u)
T2 <- f2(u)
c(mean(T1), mean(T2))  
c(sd(T1)/sqrt(m), sd(T2)/sqrt(m))

## -----------------------------------------------------------------------------
M <- 1e4
g = function (x) {
  x ^ 2 / sqrt(2*pi) * exp(-x^2/2) * (x > 1)
}

f1 = function (x) {
  dnorm(x, mean = 1.5) * (x > 1)
}
rf1 = function () {
  rnorm(M, mean = 1.5)
}
is_norm = function () {
  xs = rf1()
  return(mean(g(xs)/f1(xs), na.rm = TRUE))  
}
(theta1 = is_norm())

## -----------------------------------------------------------------------------
M <- 10000
k <- 5 
r <- M/k 
N <- 50 
T2 <- numeric(k)
est <- matrix(0, N, 2)
g<-function(t)(1-exp(-1))/(1+t^2)*(t>0)*(t<1)
for (i in 1:N) {
  est[i, 1] <- mean(g(runif(M)))
  for(j in 1:k)T2[j]<-mean(g(runif(M/k,(j-1)/k,j/k)))
  est[i, 2] <- mean(T2)
}
apply(est,2,mean)
apply(est,2,sd)

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
m <- 1000
mu1 <- 0
s1<-0
UCL <- USL <- numeric(m)
est <- matrix(0, m, 2)
for (i in 1:m) {
  y <- rnorm(n,mean = mu1,sd =2)
  est[i, 1] <- mean(y)
  est[i, 2] <- sqrt(var(y))* qt(1-alpha/2,df = n-1) / sqrt(n)
  UCL[i] <- est[i, 1] + est[i, 2]
  USL[i] <- est[i, 1] - est[i, 2]
  s1<-s1+(USL[i] < mu1 && mu1 < UCL[i] )
}
s1/m

mu2 <- 2
s2 <- 0
UCL <- USL <- numeric(m)
est <- matrix(0, m, 2)
for (i in 1:m) {
  y <- rchisq(n,df = 2)
  est[i, 1] <- mean(y)
  est[i, 2] <- sqrt(var(y))* qt(1-alpha/2,df = n-1) / sqrt(n)
  UCL[i] <- est[i, 1] + est[i, 2]
  USL[i] <- est[i, 1] - est[i, 2]
  s2<-s2+(USL[i] < mu2 && mu2 < UCL[i] )
}
s1/m
s2/m

## -----------------------------------------------------------------------------
P6A <- function(seed){
  set.seed(100)
  num<-c(50,100,200,500,1000) 
  M<-10000
  er<-NULL
  for (n in num){
    cv<-qt(0.975,n-1)
    er1 <- 0
    er2 <- 0
    er3 <- 0
    for(i in 1:M){
      x<-rchisq(n,1)
      m<-mean(x)
      se<-sqrt(var(x))
      P<- abs((m-1)*sqrt(n)/se)>=cv
      er1 <- er1+P
    }
    er1 <- er1/M
    
    for(i in 1:M){
      x<-runif(n,0,2)
      m<-mean(x)
      se<-sqrt(var(x))
      P<- abs((m-1)*sqrt(n)/se)>=cv
      er2 <- er2+P
    }
    er2 <- er2/M
    
    for(i in 1:M){
      x<-rexp(n,1)
      m<-mean(x)
      se<-sqrt(var(x))
      P<- abs((m-1)*sqrt(n)/se)>=cv
      er3 <- er3+P
    }
    er3 <- er3/M

    er<-cbind(er,c(er1,er2,er3))
  }
  colnames(er)<-num
  rownames(er)<-c("chi(1)","U(0,2)","exp(1)")
  return(er)                
}
P6A(1000)

## -----------------------------------------------------------------------------
M <- 1000      
alpha <- 0.1   

# 初始化结果变量
FWER_bonf <- 0  # Bonferroni cue FWER
FDR_bh <- 0    # B-H cue FDR
TPR_bh <- 0    # B-H cue TPR

#simulate
for (i in 1:M) {
  m <- 1000
  p <- rep(0, m)
  
  for (j in 1:m) {
    if (runif(1) <= 0.95) {
      p[j] <- runif(1)
    } else {
      p[j] <- rbeta(1, 0.1, 1)
    }
  }
  
  # Bonferroni
  p_bonf <- p.adjust(p, method = "bonferroni")
  
  # B-H
  p_bh <- p.adjust(p, method = "BH")
  
  
  reject_bonf <- p_bonf <= alpha
  reject_bh <- p_bh <= alpha
  
 
  FWER_bonf <- FWER_bonf + sum(reject_bonf)
  FDR_bh <- FDR_bh + sum(reject_bh & !reject_bonf)
  TPR_bh <- TPR_bh + sum(reject_bh)
}


FWER_bonf <- FWER_bonf / (M * m)
FDR_bh <- FDR_bh / (M * m)
TPR_bh <- TPR_bh / (M * m)

# result
result <- matrix(c(FWER_bonf, NA, NA,
                   NA, FDR_bh, TPR_bh),
                 nrow = 2, ncol = 3,
                 dimnames = list(c("Bonferroni", "B-H"),
                                 c("FWER", "FDR", "TPR")))
print(result)

## -----------------------------------------------------------------------------
lambda <- 2  
sample_sizes <- c(5, 10, 20)  
B <- 1000  
m <- 1000  

mean_bias <- matrix(0, nrow = m, ncol = length(sample_sizes))  # Initialise the matrix storing the mean deviation
std_error <- matrix(0, nrow = m, ncol = length(sample_sizes))  # Initialise the matrix for storing standard errors

for (i in 1:m) {
  for (j in 1:length(sample_sizes)) {
    n <- sample_sizes[j]  
    X <- rexp(n, rate = lambda)  
    X_bar <- mean(X)  
    
    # Resampling using the self-help method and calculating the MLE for each resample
    bootstrap_estimates <- replicate(B, 1/mean(sample(X, replace = TRUE)))
    mean_bias[i, j] <- mean(bootstrap_estimates) - lambda/(n-1)  
    std_error[i, j] <- sqrt(n-2) * mean(bootstrap_estimates) * (n/(n-1))  
  }
}

mean_bias_theoretical <- lambda/(sample_sizes - 1)  # theory deviation
std_error_theoretical <- lambda * sample_sizes / ((sample_sizes - 1) * sqrt(sample_sizes - 2))  

mean_bias_mean <- apply(mean_bias, 2, mean)  # Mean value of average deviation
std_error_mean <- apply(std_error, 2, mean)  # Mean value of standard errors

mean_bias_theoretical  # output theory deviation
std_error_theoretical  # Output Theory Standard Error
mean_bias_mean  # Mean deviation of output self-help method
std_error_mean  # Standard error of the mean of the output self-help method

## -----------------------------------------------------------------------------
library(boot)
library(bootstrap)
T_73 <- function(seed=1000){
  set.seed(seed)
  boot.t.ci <-
    function(x, B = 500, R = 100, level = .95, statistic){
      #compute the bootstrap t CI
      x <- as.matrix(x); n <- nrow(x)
      stat <- numeric(B); se <- numeric(B)
      boot.se <- function(x, R, f) {
        x <- as.matrix(x); m <- nrow(x)
        th <- replicate(R, expr = {i <- sample(1:m, size = m, replace = TRUE) 
        f(x[i, ])})
        return(sd(th))
      }
      for (b in 1:B) {
        j <- sample(1:n, size = n, replace = TRUE)
        y <- x[j, ]
        stat[b] <- statistic(y)
        se[b] <- boot.se(y, R = R, f = statistic)
      }
      stat0 <- statistic(x)
      t.stats <- (stat - stat0) / se
      se0 <- sd(stat)
      alpha <- 1 - level
      Qt <- quantile(t.stats, c(alpha/2, 1-alpha/2), type = 1)
      names(Qt) <- rev(names(Qt))
      CI <- rev(stat0 - Qt * se0)
    }
  
  dat <- cbind(law$LSAT, law$GPA)
  stat <- function(dat){mean(cor(dat[,1],dat[,2]))}
  ci <- boot.t.ci(dat, statistic = stat, B=2000, R=200)
  print(ci)
}
T_73()

## -----------------------------------------------------------------------------
library(boot)
hours = aircondit$hours
n = length(hours)
B = 200 
#MLE
mle.lambda = function (values) {
  return(length(values)/sum(values))
}

time.b = numeric(B) 
ts = numeric(B) 

time.hat = 1/mle.lambda(hours) 

# B times bootstrap sampler
for (b in 1:B) {
  i = sample(1:n, n, replace = TRUE)
  hours.b = hours[i]
  time.b[b] = 1/mle.lambda(hours.b)
  
  times.b2 = numeric(B)
  
  # Calculate parameter estimates for each bootstrap sample
  for (b2 in 1:B) {
    i2 = sample(1:n, n, replace = TRUE)
    hours.b2 = hours.b[i2]
    times.b2[b2] = 1/mle.lambda(hours.b2)
  }
  #Calculate the t-statistic for each bootstrap sample
  ts[b] = (time.b[b] - time.hat) / sd(times.b2)
}

#Calculate the standard error of the parameter estimates for the bootstrap sample
se.hat = sd(time.b)
alpha = 0.05;
q.probs = c(alpha/2, 1-alpha/2)

setCINames = function (object) {
  return(setNames(object, c(paste((alpha/2)*100, '%'), paste((1-alpha/2)*100, '%'))))
}

# Calculate confidence intervals for the standard normal distribution
q = qnorm(1-alpha/2)
ci.sn = time.hat + c(-1,1)*q*se.hat
(ci.sn = setCINames(ci.sn))

# basic boostrap.
qs.time.hat = quantile(x = time.b, p = q.probs)
ci.basic = rev(2*time.hat - qs.time.hat)
(ci.basic = setCINames (object = ci.basic))


(ci.percentile = qs.time.hat)

# Bootstrap t method
qs.t = quantile(x = ts, p = q.probs)
(ci.t = setCINames(rev(time.hat - qs.t*se.hat)))

## -----------------------------------------------------------------------------
set.seed(seed=100)
library(bootstrap)
sc<-scor
theta<-function(x){ 
  sigma<-cov(x)
  pca.sigma<-prcomp(sigma)
  theta<-pca.sigma$sdev[1]/sum(pca.sigma$sdev) 
  theta
}
n<-NROW(sc)
theta.j<- numeric(n)
for (i in 1:n){
  theta.j[i]<-theta(sc[-i,])
}
theta.hat<-theta(sc)
bias<-(n-1)*(mean(theta.j)-theta.hat) #Bias
se<-sqrt((n-1)*var(theta.j)) #Se
round(c(bias,se),3)
list(Bias.jack=bias,SE.jack=se)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- matrix(0,n,n)
for (i in 1:n) {
  
  for (j in (i+1):n-1){
    y <- magnetic[-c(i,j)]
    x <- chemical[-c(i,j)]
    
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[c(i,j)]
    e1[i,j] <- mean((magnetic[c(i,j)] - yhat1)^2)
    
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[c(i,j)] + J2$coef[3] *
      chemical[c(i,j)]^2
    e2[i,j] <- mean((magnetic[c(i,j)] - yhat2)^2)
    
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(i,j)]
    yhat3 <- exp(logyhat3)
    e3[i,j] <-mean(( magnetic[c(i,j)] - yhat3)^2)
    
    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(i,j)])
    yhat4 <- exp(logyhat4)
    e4[i,j] <- mean((magnetic[c(i,j)] - yhat4)^2)
    
  }
  
}
c(2*sum(e1)/(n*(n-1)),2*sum(e2)/(n*(n-1)),2*sum(e3)/(n*(n-1)),2*sum(e4)/(n*(n-1)))

## -----------------------------------------------------------------------------
attach(chickwts) # Attaching the chickwts dataset to the current environment

# Body weights of chickens fed with the two feed types were extracted separately and ranked
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"])) 

detach(chickwts) 

library(cramer) 

rep <- 1000 

z <- c(x, y) 
n1 <- length(x)
n2 <- length(y)
n <- n1 + n2

reps <- numeric(rep) 

T.hat <-  cramer.test(x, y)$statistic


for (i in 1:rep) {
  # A set of n1 samples is randomly selected without playback from z
  k <- sample(1:n, n1, replace = FALSE)
  z1 <- z[k]
  z2 <- z[-k]
  reps[i] <- cramer.test(z1, z2)$statistic
}

p.hat <-  mean(abs(T.hat) < abs(reps))
p.hat

hist(reps, main = "", freq = FALSE, xlab = "T (P = 0.75)", breaks = "scott")
points(T.hat, 0, cex = 1, pch = 16)

## -----------------------------------------------------------------------------
set.seed(100)
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  weight <- round(10 * (length(x)/(length(x)+length(y))))
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # Calculate the number of X's greater than Y and the number of Y's less than Y, and the number of Y's greater than X and the number of Y's less than X, based on the extreme values of the samples X and Y
  return(as.integer(min(c(outx, outy)) > weight || max(c(outx, outy)) > (10-weight) )) #return 1 reject or 0 (receive H0)
}

R <- 10000
n1 <- 10
n2 <- 20
weight <- round(10 * (n1/(n1+n2)))
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)

z <- c(x,y)
K <- 1:(n1 + n2)
n <- length(x)
jg <- numeric(R)
for (i in 1:R) {
  k <- sample(K, size = n, replace = FALSE)
  x1 <- z[k];y1 <- z[-k]
  x1 <- x1 - mean(x1) #centered by sample mean
  y1 <- y1 - mean(y1)
  jg[i] <- count5test(x1, y1)
}

q <- mean(jg) #Calculate the mean of the vector of experimental results jg, indicating the proportion of rejections of the null hypothesis
round(c(q),3) #A simulation experiment for comparing the difference between two normally distributed samples was performed by rejecting H0

## -----------------------------------------------------------------------------
ex1 <- function(N, b1, b2, b3, f_0) {
  X1 <- rbinom(N, 1, 1)
  X2 <- rexp(N, 1)
  X3 <- rbinom(N, 1, 0.5)
  #change traverse
  a <- log(f_0) - log(1 - log(f_0)) - b1 * mean(X1) - b2 * mean(X2) - b3 * mean(X3)
  
  return(a)
}

N <- 1e6
b1 <- 0
b2 <- 1
b3 <- -1
f_0_values <- c(0.1, 0.01, 0.001, 0.0001)

a_values <- c()
log_f_0_values <- c()

for (f_0 in f_0_values) {
  a <- ex1(N, b1, b2, b3, f_0)
  a_values <- c(a_values, a)
  log_f_0_values <- c(log_f_0_values, log(f_0))
}

plot(log_f_0_values, a_values, xlab = "log(f0)", ylab = "a", main = "log(f0) vs a")

## -----------------------------------------------------------------------------
Metropolis <- function(n, sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  mu <- 1
  b <- 1
  u <- runif(N,-.5,.5)
  R <- mu-b*sign(u)*log(1-2*abs(u))
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (R[i] <= (dt(y, n) / dt(x[i-1], n)))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      } }
  return(list(x=x, k=k))
}
n <- 4 #t
N <- 2000
sigma <- c(.05, .5, 2, 5)
x0 <- 25
rp1 <- Metropolis(n, sigma[1], x0, N)
rp2 <- Metropolis(n, sigma[2], x0, N)
rp3 <- Metropolis(n, sigma[3], x0, N)
rp4 <- Metropolis(n, sigma[4], x0, N)
print(c(1-rp1$k/N, 1-rp2$k/N, 1-rp3$k/N, 1-rp4$k/N))

## -----------------------------------------------------------------------------
m = 5000
burn = 1000

x = matrix(0, m, 2)

rho = 0.9
mu1 = 0
mu2 = 0
sigma1 = 1
sigma2 = 1
#Calculate the standard deviation to ensure that the generated data meets the set correlation
s1 = sqrt(1-rho^2)*sigma1
s2 = sqrt(1-rho^2)*sigma2
mean12 = function (x2) mu1 + rho*sigma1/sigma2*(x2 - mu2)
mean21 = function (x1) mu2 + rho*sigma2/sigma1*(x1 - mu1)
#Let the initial sample point be the mean
x[1,] = c(mu1, mu2)
#Generate residual samples
for (i in 2:m) {
  xt = x[i-1,]
  xt[1] = rnorm(1, mean12(xt[2]), s1)
  xt[2] = rnorm(1, mean21(xt[1]), s2)
  x[i,] = xt
}
#Discard the initial burn sample points and keep the remaining samples.
x = x[(burn+1):m,]

x = data.frame(x)
lin.reg = lm(X1 ~ X2, data = x)

par(mfrow=c(1,2))
plot(x, cex = 0.5, main = "生成数据")
hist(lin.reg$residuals, main = "线性模型残差")
par(mfrow=c(1,1))

## -----------------------------------------------------------------------------
rayleigh_pdf <- function(x, sigma) {
  if (any(x < 0)) return(0)
  stopifnot(sigma > 0)
  return(x / sigma^2 * exp(-x^2 / (2 * sigma^2)))
}

# Metropolis-Hastings
metropolis_hastings_sampler <- function(n, sigma, initial_value) {
  samples <- numeric(n)
  x <- initial_value
  
  for (i in 1:n) {
    y <- rnorm(1, mean = x, sd = 1)

    alpha <- rayleigh_pdf(y, sigma) / rayleigh_pdf(x, sigma)
    u <- runif(1)
    if (u <= alpha) {
      x <- y
    }
    samples[i] <- x
  }
  
  return(samples)
}

# Gelman-Rubin
gelman_rubin <- function(chains) {
  m <- length(chains) 
  n <- min(sapply(chains, length))
  # Intercept all chains to the same length
  chains <- lapply(chains, function(x) x[1:n])
  # Calculate the mean value of each chain
  chain_means <- sapply(chains, mean)
  # Calculate the mean of all chains
  all_means <- colMeans(do.call(rbind, chains))

  B <- n * var(chain_means)
  
  W <- mean(sapply(chains, var))
  
  var_hat <- ((n - 1) / n) * W + (1 / n) * B
  #result
  R_hat <- sqrt(var_hat / W)
  
  return(R_hat)
}


sigma <- 4  
initial_value <- 1  
target_R_hat <- 1.2  # Gelman-Rubin

chains <- list()

while (TRUE) {
  samples <- metropolis_hastings_sampler(n = 10000, sigma = sigma, initial_value = initial_value)
  
  chains <- c(chains, list(samples))
  
  # Gelman-Rubin
  R_hat <- gelman_rubin(chains)
  
  # Check that convergence conditions are met
  if (!is.na(R_hat) && R_hat < target_R_hat) {
    break
  }
  
  # Update initial sample values
  initial_value <- samples[length(samples)]
}

print(paste("Final R_hat:", R_hat))

## -----------------------------------------------------------------------------
data <- matrix(c(11, 12,
                 8, 9,
                 27, 28,
                 13, 14,
                 16, 17,
                 0, 1,
                 23, 24,
                 10, 11,
                 24, 25,
                 2, 3), ncol = 2, byrow = TRUE)

MLE <- function(lambda) {
  u <- data[, 1]
  v <- data[, 2]
  
  
  log_likelihood <- sum(log(exp(-lambda * u) - exp(-lambda * v)))
  
  return(-log_likelihood) # Minimisation using negative log-likelihood function with optim function
}

a1 <- optim(par = 1, fn = MLE, method = "Brent", lower = 0, upper = 10)

cat("直接极大似然估计结果：", a1$par, "\n")

EM <- function(data, max_iter = 10000, tol = 1e-6) {
  n <- nrow(data)
  lambda <- 1
  
  for (iter in 1:max_iter) {
    lambda_prev <- lambda
    
    # E
    z <- (exp(-lambda * data[, 1]) - exp(-lambda * data[, 2])) / (exp(-lambda * data[, 1]) - exp(-lambda * data[, 2]))
    
    # M
    lambda <- optimize(f = function(x) {
      sum(z * log(1 / x) + (1 - z) * log(-exp(-x * data[, 2]) + exp(-x * data[, 1])))
    }, interval = c(0, 10), maximum = FALSE)$minimum
    
    if (abs(lambda - lambda_prev) < tol) {
      break
    }
  }
  
  return(lambda)
}


a2 <- EM(data)


cat("EM算法估计结果：", a2, "\n")

## -----------------------------------------------------------------------------
library(boot)
# Create payment matrix A
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
               2,0,0,0,-3,-3,4,0,0,
               2,0,0,3,0,0,0,-4,-4,
               -3,0,-3,0,4,0,0,5,0,
               0,3,0,-4,0,-4,0,5,0,
               0,3,0,0,4,0,-5,0,-5,
               -4,-4,0,0,0,5,0,0,6,
               0,0,4,-5,-5,0,0,0,6,
               0,0,4,0,0,5,-6,-6,0), 9, 9)
B <- A + 2

solve.game <- function(A){
  # Solving a zero-sum game with two players using the simplex method
  # Optimise player 1 first, optimise player 2 later
  
   
  min.A <- min(A)
  A <- A - min.A # v >= 0
  max.A <- max(A)
  A <- A / max(A)
  
  m <- nrow(A)
  n <- ncol(A)
  it <- n^3
  
  # Player 1 optimisation
  a <- c(rep(0, m), 1)  
  A1 <- -cbind(t(A), rep(-1, n)) 
  b1 <- rep(0, n) 
  A3 <- t(as.matrix(c(rep(1, m), 0))) 
  b3 <- 1
  sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3, maxi=TRUE, n.iter=it)
  
  # Player 2 optimisation
  a <- c(rep(0, n), 1)  
  A1 <- cbind(A, rep(-1, m))  
  b1 <- rep(0, m)  
  A3 <- t(as.matrix(c(rep(1, n), 0)))  
  b3 <- 1
  sy <- simplex(a = a, A1 = A1, b1 = b1, A3 = A3, b3 = b3, maxi = FALSE, n.iter = it)
  
  # 构造结果
  soln <- list(
    "A" = A * max.A + min.A,  # 恢复原始范围的支付矩阵A
    "x" = sx$soln[1:m],  # Player 1's strategy vector x
    "y" = sy$soln[1:n],  # Player 2's strategy vector x
    "v" = sx$soln[m + 1] * max.A + min.A,  # The value of the game v
    "optimal_point" = cbind(sx$soln[1:m], sy$soln[1:n])  
  )
  
  return(soln)
}

soln_A <- solve.game(A)
soln_B <- solve.game(B)
round(soln_A$optimal_point, 7)
round(soln_A$v, 7)
round(soln_B$optimal_point, 7)
round(soln_B$v, 7)

## -----------------------------------------------------------------------------
# In R, the unlist() function is used to convert lists into atomic vectors by combining all the elements of a list into a single vector, while the as.vector() function is mainly used to force conversions and cannot directly convert lists into atomic vectors.
#The main reason for using unlist() instead of as.vector() is that unlist() is specifically designed to flatten a list structure and concatenate its elements into a single vector. It recursively traverses the list and extracts each element, preserving their order. This feature of unlist() is critical when there are nested lists with multiple levels of nesting and you want to collapse them into a single vector.
#as.vector() is primarily used for forced conversions, where it attempts to convert objects into vectors of the appropriate pattern. It is useful for converting other data types (such as matrices or arrays) to vector form. However, when it comes to lists, as.vector() does not automatically flatten the list structure. It can only force a list to be converted to a vector of type "list", where each element of the resulting vector is a list.

## -----------------------------------------------------------------------------
#dim() function is used to retrieve or set the dimension of an object, such as a matrix or array. When applied to a vector, the dim() function returns NULL.
# In R, a vector is considered a one-dimensional object that represents a sequence of elements with the same data type. Since vectors do not have multiple dimensions like matrices or arrays, the dim() function does not have any meaningful dimensions to return.
dim(c(1, 2, 3))

## -----------------------------------------------------------------------------
# Again, TRUE, the array is a multidimensional matrix.
x <- matrix(1:4, nrow = 2)
is.matrix(x)
is.array(x)

## -----------------------------------------------------------------------------
# Forced type conversion occurs and the resulting matrix will be a character matrix
df <- data.frame(A = c(1, 2, 3), B = c("a", "b", "c"), C = c(TRUE, FALSE, TRUE))

matrix <- as.matrix(df)

print(matrix)
print(class(matrix))

## -----------------------------------------------------------------------------
#yes
#The number of rows and columns is 0
empty_df <- data.frame()

num_rows <- nrow(empty_df)
num_cols <- ncol(empty_df)

print(num_rows)  
print(num_cols) 

#rows are 0
df <- data.frame(1, seq(10))[FALSE, ]
dim(df)
#cols are 0
df <- data.frame(1, seq(10))[, FALSE]
dim(df)

## -----------------------------------------------------------------------------
#scale01 <- function(x) {
#  rng <- range(x, na.rm = TRUE)
#  (x - rng[1]) / (rng[2] - rng[1])
#}

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}
lapply(mtcars, scale01) 

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}
#As an example, the Iris dataset
iris[sapply(iris, is.numeric)]#Selection series
lapply(iris[sapply(iris, is.numeric)], scale01)

## -----------------------------------------------------------------------------
#a
df<-data.frame(x=c(1,2,3),y=c(4,4,5))
df
vapply(df,sd,FUN.VALUE = 1)

#b
df<-data.frame(x=c(1,2,3),y=c(4,4,5))
df$z<-c("A","B","C")
df
vapply(df[vapply(df,is.numeric,TRUE)],sd,1)

## ----eval=FALSE---------------------------------------------------------------
#  #R code:
#  GibbsR <- function(N,thin,a,b,n){
#    mat <- matrix(nrow = N, ncol = 2)
#    x <-rbinom(1,prob=0.5,size=n)
#    y <-rbeta(1,x+a,n-x+b)
#    for (i in 1:N) {
#      for (j in 1:thin) {
#        x <- rbinom(1,prob = y, size = n)
#        y <- rbeta(1,x+ a, n - x+ b)
#      }
#      mat[i, ] <- c(x, y)
#    }
#    mat
#  }
#  
#  #The C++ code:
#  library(Rcpp)
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  
#  // [[Rcpp::export]]
#  NumericMatrix GibbsC(int N, int thin, int a, int b,int n) {
#     NumericMatrix mat(N, 2);
#     double x = rbinom(1,0.5,n)[0];
#     double y = rbeta(1,x+a,n-x+b)[0];
#     for(int i = 0; i < N; i++) {
#       for(int j = 0; j < thin; j++) {
#         x = rbinom(1,y,n)[0];
#         y = rbeta(1,x+a,n-x+b)[0];
#       }
#       mat(i, 0) = x;
#       mat(i, 1) = y;
#     }
#     return(mat);
#  }
#  
#  library(microbenchmark)
#  tm1 <- microbenchmark::microbenchmark(
#  rnR = GibbsR(1e3,10,2,3,10),
#  rnC = GibbsC(1e3,10,2,3,10)
#   )
#   print(summary(tm1)[,c(1,3,5,6)])


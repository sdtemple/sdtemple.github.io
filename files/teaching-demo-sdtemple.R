# UW Biostatistics Seminar
# Seth Temple, sdtemple@uw.edu
# January 18, 2024



# Color blind color palettes ----------------------------------------------

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pie(rep(1, 8), col = colorBlindBlack8)

col1=colorBlindBlack8[1]
col2=colorBlindBlack8[2]
col3=colorBlindBlack8[3]
col4=colorBlindBlack8[4]

# Sodium metabolite as normally distributed -------------------------------

# help functions
?rnorm
?hist
?curve
?abline

mu = 140
sig = 2.5
x = rnorm(5000, mean=mu, sd=2.5)

# plotting
hist(x,
     breaks=50,
     xlab='Sodium (mEq per liter)',
     main=NA,
     probability=T,
     xlim=c(130,150)
     )
curve(dnorm(x,mu,sig),
      from=130,
      to=150,
      add=T,
      lwd=3
      )
# abline(v=mu,
#        col=col3,
#        lwd=4,
#        lty='dashed'
#        )
abline(v=135,
       col=col3,
       lwd=2,
       lty='dashed'
       )
abline(v=145,
       col=col3,
       lwd=2,
       lty='dashed'
)
legend("topright",
       legend=c("Adverse Risks"),
       col=c(col3),
       lty=c("dashed"),
       lwd=c(2)
)

# p-value calculation
n = 30
mud = 139
xd = rnorm(n, mud, sig)
zd = sqrt(n) * (mean(xd) - mu) / sig
pnorm(zd,0,1)


# Type Error Simulations --------------------------------------------------

# type 1 error

ctr <- 0
K <- 1000
n = 30
mu = 140
sig = 2.5
for(j in 1:K){
  x = rnorm(n,mu,sd=sig)
  z = sqrt(n) * (mean(x) - mu) / sig
  if(abs(z)>1.96){
    ctr <- ctr + 1
  }
} 
ctr/K

# alternative model diagram

?abline

vert = mu - 2 * sig
x = rnorm(n,mu,sd=sig)
y = rnorm(n,138,sig)

curve(dnorm(x,mu,sig),
      from=130,
      to=150,
      lwd=3,
      xlab='Sodium (mEq per liter)',
      y='Density',
      col=col1
)
curve(dnorm(x,138,sig),
      from=130,
      to=150,
      lwd=3,
      xlab='Sodium (mEq per liter)',
      y='Density',
      col=col2,
      add=T
)

abline(v=135,
       col=col3,
       lwd=3,
       lty='dashed'
)
legend("topright",
       legend=c("Null","Alt","Bad"),
       col=c(col1,col2,col3),
       lty=c("solid","solid","dashed"),
       lwd=c(3,3,3)
       )


# type 2 error (different models)

{
ctr <- 0
K <- 2000
n = 40
mu = 140
mua = 139. # change this and re-run
sig = 2.5
for(j in 1:K){
  x = rnorm(n,mua,sd=sig)
  z = sqrt(n) * (mean(x) - mu) / sig
  if(abs(z)>1.96){
    ctr <- ctr + 1
  }
}
ctr/K
(1-ctr/K)
}

# type 2 error (different sample sizes)

{
  ctr <- 0
  K <- 2000
  n = 40 # change this and re-run
  mu = 140
  mua = 139. 
  sig = 2.5
  for(j in 1:K){
    x = rnorm(n,mua,sd=sig)
    z = sqrt(n) * (mean(x) - mu) / sig
    if(abs(z)>1.96){
      ctr <- ctr + 1
    }
  }
  ctr/K
  (1-ctr/K)
}

# type 2 error (different level tests)

{
K <- 2000
n = 30 
mu = 140
mua = 139.
sig = 2.5
cutoff = qnorm(0.9995) # change this and re-run
ctr <- 0
for(j in 1:K){
  x = rnorm(n,mua,sd=sig)
  z = sqrt(n) * (mean(x) - mu) / sig
  if(abs(z)>cutoff){
    ctr <- ctr + 1
  }
} 
ctr/K
(1-ctr/K)
}

# # hidden slide
# 
# mu <- 140
# sig <- 2.5
# n <- 8
# x <- rnorm(1,mu,sig/sqrt(n))
# curve(dnorm(x,mu,sig/sqrt(n)),
#       from=130,
#       to=150,
#       lwd=2,
#       xlab='Sodium (mEq per liter)',
#       y='Density',
#       main="Sample distribution of average"
# )

# n <- 4
# x <- rnorm(1,mu,sig/sqrt(n))
# curve(dnorm(x,mu,sig/sqrt(n)),
#       from=130,
#       to=150,
#       lwd=2,
#       add=T,
#       col=col2,
#       lty="dotted"
# )
# 
# n <- 1
# x <- rnorm(1,mu,sig/sqrt(n))
# curve(dnorm(x,mu,sig/sqrt(n)),
#       from=130,
#       to=150,
#       lwd=2,
#       add=T,
#       col=col3,
#       lty="dashed"
# )
# 
# abline(v=139,lty="dashed")
# 
# legend("topright",
#        legend=c("8","4","1"),
#        col=c(col1,col2,col3),
#        lty=c("solid","dotted","dashed"),
#        lwd=c(2,2,2),
#        title="Sample Size"
# )



# Genome-wide associations ------------------------------------------------

?lm
?summary
?rmultinom

# additive genetic effect
n = 50
eps = rnorm(n,sd=sig)
x = rmultinom(n, 2, rep(1/3,2))
x = x[1,]
beta1 = 3
beta0 = 137
y = beta0 + beta1 * x + eps

plot(x,
     y,
     pch=20,
     xlab="Genotype",
     ylab="Sodium Level",
     ylim=c(130,150),
     main="Additive genetic effect"
     )
abline(a=beta0,
       b=beta1,
       lwd=2,
       col=col2
       )
summary(lm(y~x))

# no genetic effect
n = 50
eps = rnorm(n,sd=sig)
x = rmultinom(n, 2,rep(1/3,2))
x = x[1,]
beta1 = 0
beta0 = 140
y = beta0 + beta1 * x + eps

plot(x,
     y,
     pch=20,
     xlab="Genotype",
     ylab="Sodium Level",
     ylim=c(130,150),
     main="No genetic effect"
)
abline(a=beta0,
       b=beta1,
       lwd=2,
       col=col2
)

summary(lm(y~x))

# multiple testing

pvals <- c()
n = 50
beta1 = 0
beta0 = 140
sig=2.5
ctr=0
K <- 1000
for(k in 1:K){
  eps = rnorm(n,sd=sig)
  x = rmultinom(n, 2,rep(1/3,2))
  x = x[1,]
  y = beta0 + beta1 * x + eps
  currpval = summary(lm(y~x))$coefficients[,4][2]
  if(currpval<0.01){
    ctr <- ctr + 1
  }
  pvals <- c(pvals,currpval)
}

plot(1:K,
     -log10(pvals),
     pch=20,
     xlab="Genetic Marker",
     ylab="-log10 p-value",
     main="Is sodium metabolite associated w/ genotypes?",
     ylim=c(0,7)
     )
abline(h=2,col=col2,lwd=2,lty="dashed")
abline(h=5,col=col3,lwd=2,lty="dotted") # comment out for 1 slide

legend("topright",
       legend=c("0.01","0.01 / 1000"),
       col=c(col2,col3),
       lty=c("dashed","dotted"),
       lwd=c(2,2),
       title="p-value cutoff"
)

legend("topright",
       legend=c("0.01"),
       col=c(col2),
       lty=c("dashed"),
       lwd=c(2),
       title="p-value cutoff"
)

# Computing Z test statistic ----------------------------------------------

mu = 140
sig = 2.5
x = rnorm(5000, mean=mu, sd=2.5)

# plotting
hist(x,
     breaks=50,
     xlab='Sodium (mEq per liter)',
     main=NA,
     probability=T,
     xlim=c(130,150)
)
curve(dnorm(x,mu,sig),
      from=130,
      to=150,
      add=T,
      lwd=3
)

# plotting
hist(x-140,
     breaks=50,
     xlab='Sodium (mEq per liter)',
     main=NA,
     probability=T,
     xlim=c(-10,10)
)
curve(dnorm(x,0,sig),
      from=-10,
      to=10,
      add=T,
      lwd=3
)

# plotting
hist((x-140)/2.5,
     breaks=50,
     xlab='Sodium (mEq per liter)',
     main=NA,
     probability=T,
     xlim=c(-3,3)
)
curve(dnorm(x,0,1),
      from=-3,
      to=3,
      add=T,
      lwd=3
)

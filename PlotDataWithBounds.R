####### R code for Plotting Data with Bounds #######
####### Version 1.0: (8/30/2017)     #######
####### By Alan Aw                   #######

#### Helper Functions ####
little.a <- function(n) {
  harm.sum = 0
  for (i in 1:(n-1)) {
    harm.sum = harm.sum + 1/i
  }
  return(harm.sum)
}

big.a <- function(n) {
  harm.sum.square = 0
  for (i in 1:(n-1)) {
    harm.sum.square = harm.sum.square + 1/(i^2)
  }
  return(harm.sum.square)
}

### For Tajima's D ### (see Simonsen, Churchill, Aquadro, 1995)
v.T <- function(n) {
  (2*(n^2+n+3)/(9*n*(n-1))-(n+2)/(little.a(n)*n) + big.a(n)/(little.a(n))^2) / (little.a(n)^2+big.a(n))
}

u.T <- function(n) {
  (((n+1)/(3*(n-1))-1/little.a(n))/little.a(n)) - v.T(n) 
} ###

#### Tajima's D ####

### Upper Bound ###
f.TajimaD.max <- function(n,S,j,k) {
  return((f.thetaPi.max(n,S,j,k) - S/little.a(n)) / sqrt(u.T(n)*S + v.T(n)*S*S))
}

### Lower Bound ###
f.TajimaD.min <- function(n,S,j,k) {
  (f.thetaPi.min(n,S,j,k) - S/little.a(n)) / sqrt(u.T(n)*S + v.T(n)*S*S)
}

## Let's plot dependence of
## Tajima's Pi on k ##
vec.k = c(0:S) 
vec.thetaPi.max = f.thetaPi.max(n,S,j,vec.k)
plot(vec.k/S,vec.thetaPi.max, xlab="Frequency of Derived Singletons", ylab="",col="blue",ylim=c(0,50),main=paste("n =",n,", S =", S),type="n")
lines(vec.k/S, vec.thetaPi.max, xlim=range(vec.k), ylim=range(vec.thetaPi.max), col="blue4",lwd=2,pch=16)
vec.thetaPi.min = f.thetaPi.min(n,S,j,vec.k) 
points(vec.k/S,vec.thetaPi.min,col="red",type="n") 
lines(vec.k/S, vec.thetaPi.min, xlim=range(vec.k), ylim=range(vec.thetaPi.min), col="darkred",lwd=2,pch=16)

# Copy and paste: expression("Tajima's"~pi), Frequency of Derived Singletons

## Let's plot dependence of
## Tajima's D on frequency of singletons, k ##
vec.k = c(0:S)
vec.TajimaD.max = f.TajimaD.max(n,S,j,vec.k)
plot(vec.k/S,vec.TajimaD.max, xlab="", ylab="",col="blue",ylim=c(-4,8),main=paste("n =",n,", S =",S),type="n")
lines(vec.k/S, vec.TajimaD.max, xlim=range(vec.k), ylim=range(vec.TajimaD.max), col="blue4",lwd=2,pch=16)
vec.TajimaD.min = f.TajimaD.min(n,S,j,vec.k)   
points(vec.k/S,vec.TajimaD.min,col="red",type="n") 
lines(vec.k/S, vec.TajimaD.min, xlim=range(vec.k), ylim=range(vec.TajimaD.min), col="darkred",lwd=2,pch=16)

## Let's plot dependence of
## Fu and Li's F on frequency of singletons, k ##
vec.k = c(0:S)
vec.FuLiF.max = f.FuLiF.max(n,S,j,vec.k)
plot(vec.k/S,vec.FuLiF.max, xlab="Frequency of Derived Singletons", ylab="",col="blue",ylim=c(-12,6),main=paste("n =",n,", S =",S),type="n")
lines(vec.k/S, vec.FuLiF.max, xlim=range(vec.k/S), ylim=range(vec.FuLiF.max), col="blue4",lwd=2,pch=16)
vec.FuLiF.min = f.FuLiF.min(n,S,j,vec.k)   
points(vec.k/S,vec.FuLiF.min,col="red",type="n") 
lines(vec.k/S, vec.FuLiF.min, xlim=range(vec.k/S), ylim=range(vec.FuLiF.min), col="darkred",lwd=2,pch=16)

# copy and paste: Fu and Li's F, "Frequency of Derived Singletons"

# Paste this in: "Frequency of Derived Singletons" "Tajima's D"
## Let's add the simulated data ##
points(k.vec,diffStat.vec,col="black",pch="o")
#====================
df = read.table("sfs.txt",header=FALSE)
#=====================
df = read.table("fuliF-bottleneckN500S100.txt",header=FALSE)
points(as.vector(df$V1[-1]),as.vector(df$V2[-1]),col="black",pch="o")

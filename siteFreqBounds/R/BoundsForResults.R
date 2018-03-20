####### R code for Theta Bounds and Difference Statistics   #######
####### Version 3.0: (3/17/2018)                            #######
####### By Alan Aw, Rayne Hernandez, Younes Bensouda Mourri #######

##### Bounds are mathematical and depend only on data, so 
##### works for any parametric molecular evolutionary model 
##### (infinite sites, infinite alleles) #####

#================================================================================================

#' Helper Function computes truncated harmonic sum
#' 
#' Comment
#' 
#' @param n where to truncate the sum
#' @return the truncated harmonic sum
#' @export
#' @examples

little.a <- function(n) {
  harm.sum = 0
  for (i in 1:(n-1)) {
    harm.sum = harm.sum + 1/i
  }
  return(harm.sum)
}

#' Helper Function computes truncated sum of inverse squares
#' 
#' Comment
#' 
#' @param n where to truncate the sum
#' @return the truncated sum
#' @export
#' 

big.a <- function(n) {
  harm.sum.square = 0
  for (i in 1:(n-1)) {
    harm.sum.square = harm.sum.square + 1/(i^2)
  }
  return(harm.sum.square)
}

## For Tajima's D ## (see Simonsen, Churchill, Aquadro, 1995)
#' Helper Function 
#' 
#' Comment
#' 
#' @param n where to truncate the sum
#' @return the function
#' @export
#' 

v.T <- function(n) {
  (2*(n^2+n+3)/(9*n*(n-1))-(n+2)/(little.a(n)*n) + big.a(n)/(little.a(n))^2) / (little.a(n)^2+big.a(n))
}

## For Tajima's D ## (see Simonsen, Churchill, Aquadro, 1995)
#' Helper Function 
#' 
#' Comment
#' 
#' @param n where to truncate the sum
#' @return the function
#' @export
#' 

u.T <- function(n) {
  (((n+1)/(3*(n-1))-1/little.a(n))/little.a(n)) - v.T(n) 
}

### For Fu and Li's F ### (see Simonsen, Churchill, Aquadro, 1995)
#' Helper Function 
#' 
#' Comment
#' 
#' @param n where to truncate the sum
#' @return the function
#' @export
v.F <- function(n) {
  ((2*n^3 + 110*n^2 -255*n + 153)/(9*n^2*(n-1))+(2*(n-1)*little.a(n))/(n^2)-(8*big.a(n))/n) / (little.a(n)^2+big.a(n))
}

#' Helper Function 
#' 
#' Comment
#' 
#' @param n where to truncate the sum
#' @return the function
#' @export

u.F <- function(n) {
  ((4*n^2+19*n+3-12*(n+1)*little.a(n+1))/(3*n*(n-1)))/little.a(n) - v.F(n)
}

### Pre-simulation step ###
## HOW MANY REPS DO YOU WANT TO SIMULATE? ##
#R = 1000
#n = 500 # equals to first numeric argument in system formula
#theta = 5.0 # equals to "-t" in system formula
#S = 100 # equals to "-s" in system formula

## (Optional) LET'S CREATE A DATA FRAME TO STORE SIMULATION OUTPUT! ##  
# df = data.frame(matrix(ncol=2,nrow=R,dimnames=list(c(),c("k","diffStat"))),stringsAsFactors=F)

#' Helper Function 
#' 
#' Comment
#' 
#' @param n comment
#' @param S comment
#' @param theta comment
#' @param scenario comment
#' @param ... comment
#' @return the function
#' 
#' @export

getCommand <- function(n, S, theta, scenario = 'neutral', ...) {
  print(getwd())
  if (scenario == 'neutral') {
    # result = paste("cd 00_pkg_src/siteFreqBounds/src/msdir && ", "./ms", n, "1 -t", theta, "-s", S, "> sfs.txt", "&& cp sfs.txt ../../../", sep = " ")
    result = paste("cd siteFreqBounds/src/msdir && ", "./ms", n, "1 -t", theta, "-s", S, "> sfs.txt", "&& cp sfs.txt ../../../", sep = " ")
  } else if (scenario == 'expansion') {
    # result = paste("cd 00_pkg_src/siteFreqBounds/src/msdir && ", "./ms", n, "1 -t", theta, "-s", S, "-G 6.93 -eG 0.2 0.0", "> sfs.txt", "&& cp sfs.txt ../../../", sep = " ")
    result = paste("cd siteFreqBounds/src/msdir && ", "./ms", n, "1 -t", theta, "-s", S, "-G 6.93 -eG 0.2 0.0", "> sfs.txt", "&& cp sfs.txt ../../../", sep = " ")
  } else if (scenario == 'bottleneck') {
    # result = paste("cd 00_pkg_src/siteFreqBounds/src/msdir && ", "./ms", n, "1 -t", theta, "-s", S, "-eN 0.2 5.0", "> sfs.txt", "&& cp sfs.txt ../../../", sep = " ")
    result = paste("cd siteFreqBounds/src/msdir && ", "./ms", n, "1 -t", theta, "-s", S, "-eN 0.2 5.0", "> sfs.txt", "&& cp sfs.txt ../../../", sep = " ")
  }
  # print(result)
  result
}

#' Helper Function 
#' 
#' Comment
#' 
#' @param R comment
#' @param n comment
#' @param S comment
#' @param theta comment
#' @param f comment
#' @param scenario comment
#' @param ... comment
#' @return the function
#' 
#' @export

runCommand <- function(R, n, S, theta, f = 'FuLi', scenario = 'neutral',...) {
  k.vec = c()
  diffStat.vec = c()
  command = getCommand(n, S, theta, scenario)
  for (r in 1:R) { 
    system(command)
    system("sed '1,7d' sfs.txt > sfs-modified.txt")
    system("sed -e 's/./,&/g' sfs-modified.txt > sfs.csv")
    sfs.matrix = read.csv("sfs.csv",header=FALSE)
    sfs.matrix = sfs.matrix[-1] 
    sum(sfs.matrix$V2)
    alleleFreq.vector=colSums(sfs.matrix)
    
    ### This step produces the SFS vector ###
    sfs.vec = c()
    for (i in 1:(n-1)) {
      xi.i = 0
      for (j in 1:S) {
        if (alleleFreq.vector[j] == i) xi.i = xi.i+1
      }
      sfs.vec = c(sfs.vec,xi.i)
    }
    
    ### This step computes Tajima's Pi / Tajima's D / Fu and Li's F and k for the vector ###
    kOverS = sfs.vec[1] / S # get frequency of derived singletons
    iOmega.vec = c(1:(n-1))
    for (i in 1:(n-1)) {
      if(f == 'FuLi'){
        iOmega.vec[i] = ifelse(i==1,1,0)*(2/n-1) + ifelse(i==1,0,1)*((n-i)*i/choose(n,2)) # Fu and Li's F (TOGGLE; uncomment to select)
      } else if (f == 'TajD'){
        iOmega.vec[i] = (n-i)*i/choose(n,2) - 1/(little.a(n)) # Tajima's D (TOGGLE; uncomment to select)
      } else if (f == 'TajPi') {
        iOmega.vec[i] = i*(n-i) / choose(n,2)
      }
      
    }
    if(f == 'FuLi'){
      fuliF=sum(iOmega.vec*sfs.vec)/sqrt(u.F(n)*S + v.F(n)*S*S) # Fu and Li's F (TOGGLE; uncomment to select)
    } else if (f == 'TajD'){
      tajimaD=sum(iOmega.vec*sfs.vec)/sqrt(u.T(n)*S + v.T(n)*S*S) # Tajima's D (TOGGLE; uncomment to select)
    } else if (f == 'TajPi') {
      tajimaPi=sum(iOmega.vec*sfs.vec) # Tajima's Pi (TOGGLE; uncomment to select)
    }
    ### This step adds simulated data to output list ###
    k.vec = c(k.vec,kOverS)
    if( f == 'FuLi'){
      diffStat.vec=c(diffStat.vec,fuliF) # Fu and Li's F (TOGGLE; uncomment to select)
    } else if (f == 'TajD'){
      diffStat.vec=c(diffStat.vec,tajimaD) # Tajima's D (TOGGLE; uncomment to select)
    } else if (f == 'TajPi'){
      diffStat.vec=c(diffStat.vec,tajimaPi) # Tajima's Pi (TOGGLE; uncomment to select)
    }
  }
  simulation.output = data.frame(k.vec,diffStat.vec)
  write.table(simulation.output, file = "output.txt", row.names = FALSE)
}  


#================================================================================================
#### Paramaters ####
#S = 100 # Segregating Sites
#n = 500 # Individuals
#j = 1 # Component of SFS, must be in [1,n-1]
#k = 5 # Value of \xi_j, must be in [0,S]
####

###########################################
######### Theta Estimator Bounds ##########
###########################################

#### Tajima's Pi ####
#' Helper Function 
#' 
#' Comment
#' 
#' @param n comment
#' @param j comment
#'
#' @return the function
#' 
#' @export

f.max <- function(n,j) { # Helper function: get upper bound
  f.max.out = 0 
  for (i in 1:(n-1)) {
    if (i!=j & i*(n-i) > f.max.out) {
      f.max.out = i * (n-i)
    }
  }
  return(f.max.out)
}

#' Helper Function 
#' 
#' Comment
#' 
#' @param n comment
#' @param j comment
#' 
#' @return the function
#' 
#' @export

f.min <- function(n,j) { # Helper function: get lower bound
  f.min.out = n * n 
  for (i in 1:(n-1)) {
    if (i!=j & i*(n-i) < f.min.out) {
      f.min.out = i * (n-i)
    }
  }
  return(f.min.out)
}

### Upper Bound ###
#' Helper Function 
#' 
#' Comment
#' 
#' @param n comment
#' @param S comment
#' @param j comment
#' @param k comment
#' @return the function
#' @export 

f.thetaPi.max <- function(n,S,j,k) { 
  (k*j*(n-j) + (S-k)*f.max(n,j)) / choose(n,2)
}

### Lower Bound ###
#' Helper Function 
#' 
#' Comment 
#' 
#' @param n comment
#' @param S comment
#' @param j comment
#' @param k comment
#' @return the function
#' @export

f.thetaPi.min <- function(n,S,j,k) {
  (k*j*(n-j) + (S-k)*f.min(n,j)) / choose(n,2)
} 

### Example Plots ###
## Initialize some parameters ##
#R=100
#n=500
#S=100
#j=1
#k=5

## Let's plot dependence of
## Tajima's Pi on k ##
#' Helper Function 
#' 
#' Comment 
#' @param R comment
#' @param n comment
#' @param S comment
#' @param j comment
#' @param theta comment
#' @param scenario comment
#' @param plot_points comment 
#' @param ... comment
#' @return the function
#' 
#' @export

plot_tajima_Pi_k <- function(R,n,S,j,theta,scenario,plot_points=FALSE, ...){
  vec.k = c(0:S) 
  vec.thetaPi.max = f.thetaPi.max(n,S,j,vec.k)
  vec.thetaPi.min = f.thetaPi.min(n,S,j,vec.k) 
  plot(vec.k/S,vec.thetaPi.max, xlab="Frequency of Derived Singletons", ylab=expression("Tajima's"~pi),col="blue",
       type='l',lwd=2,ylim=range(c(vec.thetaPi.min, vec.thetaPi.max)),main=paste("n =",n,", S =", S),pch="o")
  #lines(vec.k/S, vec.thetaPi.max, xlim=range(vec.k), ylim=range(vec.thetaPi.max), col="blue4",lwd=2,pch=16)
  #points(vec.k/S,vec.thetaPi.min,col="red",pch="o") 
  lines(vec.k/S, vec.thetaPi.min, xlim=range(vec.k), ylim=range(vec.thetaPi.min), col="red",lwd=2,pch=16)
  if(plot_points == TRUE){
    runCommand(R, n, S,theta, f='TajPi', scenario)
    df = read.table("output.txt",header=FALSE)
    points(as.vector(df$V1[-1]),as.vector(df$V2[-1]),col="black",pch="o")
  }
}

# plot_tajima_Pi_k(R,n,S,j,theta,scenario = "bottleneck",plot_points = TRUE)
# =====================================================================================================

## Let's plot dependence of
## Tajima's Pi on n ##
# n.lower=3
# n.upper=50

#' Helper Function 
#' 
#' Comment
#' 
#' @param n.lower comment
#' @param n.upper comment
#' @param S comment
#' @param j comment
#' @param k comment
#' @return the function
#' @export

plot_tajima_Pi_n <- function(n.lower, n.upper,S,j, k){
  vec.n = c(n.lower:n.upper)
  vec.thetaPi.max=c(1:length(vec.n))
  vec.thetaPi.min=c(1:length(vec.n))
  for (i in 1:length(vec.n)) {
    vec.thetaPi.max[i]=f.thetaPi.max(i+n.lower-1,S,j,k)
    vec.thetaPi.min[i]=f.thetaPi.min(i+n.lower-1,S,j,k) 
  }
  plot(vec.n,vec.thetaPi.max, xlab="No. of Sampled Individuals", ylab=expression("Tajima's"~pi),col="purple",ylim=range(c(vec.thetaPi.min, vec.thetaPi.max)),main=paste("k =",k,", S =",S),pch="o")
  lines(vec.n, vec.thetaPi.max, xlim=range(vec.n), ylim=range(vec.thetaPi.max), col="darkslateblue",lwd=2,pch=16)
  points(vec.n,vec.thetaPi.min,col="green",pch="o") 
  lines(vec.n, vec.thetaPi.min, xlim=range(vec.n), ylim=range(vec.thetaPi.min), col="darkgreen",lwd=2,pch=16)
}

# plot_tajima_Pi_n(n.lower,n.upper,S,j, k)
# =====================================================================================================

## Let's plot dependence of
## Tajima's Pi on S ##
#' Helper Function 
#' 
#' Comment
#' 
#' @param n comment
#' @param j comment
#' @param k comment
#' @return the function
#' @export

plot_tajima_Pi_S <- function(n,j,k){
  vec.S = seq(from=k,to=100,by=1) 
  vec.thetaPi.max = f.thetaPi.max(n,vec.S,j,k)
  vec.thetaPi.min = f.thetaPi.min(n,vec.S,j,k)
  plot(vec.S,vec.thetaPi.max, xlab="No. of Segregating Sites", ylab=expression("Tajima's"~pi),col="orange",ylim=range(c(vec.thetaPi.min, vec.thetaPi.max)),main=paste("k =",k,", n =",n),pch="o")
  lines(vec.S, vec.thetaPi.max, xlim=range(vec.S), ylim=range(vec.thetaPi.max), col="orange4",lwd=2,pch=16)
  points(vec.S,vec.thetaPi.min,col="pink",pch="o") 
  lines(vec.S, vec.thetaPi.min, xlim=range(vec.S), ylim=range(vec.thetaPi.min), col="pink4",lwd=2,pch=16)
}

# plot_tajima_Pi_S(n, S, j,k)

# =====================================================================================================

###########################################
########## Difference Statistics ##########
###########################################

#### Helper Functions ####



#### Tajima's D ####

### Upper Bound ###
#' Helper Function 
#' 
#' Comment
#' 
#' @param n comment
#' @param S comment
#' @param j comment
#' @param k comment
#' @return the function
#' @export
#' 
f.TajimaD.max <- function(n,S,j,k) {
  return((f.thetaPi.max(n,S,j,k) - S/little.a(n)) / sqrt(u.T(n)*S + v.T(n)*S*S))
}

### Lower Bound ###
#' Helper Function 
#' 
#' Comment
#' 
#' @param n comment
#' @param S comment
#' @param j comment
#' @param k comment
#' @return the function
#' @export
#' 
f.TajimaD.min <- function(n,S,j,k) {
  (f.thetaPi.min(n,S,j,k) - S/little.a(n)) / sqrt(u.T(n)*S + v.T(n)*S*S)
}

## Let's plot dependence of
## Tajima's D on frequency of singletons, k ##
#' Helper Function 
#' 
#' Comment 
#' 
#' @param R comment
#' @param n comment
#' @param S comment
#' @param j comment
#' @param theta comment
#' @param scenario comment
#' @param plot_points comment
#' @param ... comment
#' @return the function
#' @export
#' 
plot_tajima_D_k <- function(R,n,S,j,theta,scenario,plot_points = FALSE,...){
  vec.k = c(0:S)
  vec.TajimaD.max = f.TajimaD.max(n,S,j,vec.k)
  vec.TajimaD.min = f.TajimaD.min(n,S,j,vec.k)   
  plot(vec.k/S,vec.TajimaD.max, xlab="Frequency of Derived Singletons", ylab="Tajima's D",col="blue",
       type='l',lwd=2,ylim=range(c(vec.TajimaD.max, vec.TajimaD.min)),main=paste("n =",n,", S =",S),pch="o")
  #lines(vec.k/S, vec.TajimaD.max, xlim=range(vec.k), ylim=range(vec.TajimaD.max), col="blue4",lwd=2,pch=16)
  #points(vec.k/S,vec.TajimaD.min,col="red",pch="o") 
  lines(vec.k/S, vec.TajimaD.min, xlim=range(vec.k), ylim=range(vec.TajimaD.min), col="red",lwd=2,pch=16)
  if(plot_points == TRUE){
    runCommand(R,n,S,theta,f = 'TajD',scenario)
    df = read.table("output.txt",header=FALSE)
    points(as.vector(df$V1[-1]),as.vector(df$V2[-1]),col="black",pch="o")
  }
}
# plot_tajima_D_k(R,n,S,j,theta,scenario,TRUE)

#################====================

#n = 50 
#S = 10 
#j = 1
#plot_tajima_D_k(n,S,j)

## Let's plot dependence of
## Tajima's D on S ##
#' Helper Function 
#' 
#' Comment
#' 
#' @param n comment
#' @param j comment
#' @param k comment
#' @return the function
#' @export
#' 
plot_tajima_D_S <- function(n,j,k){
  vec.S = seq(from=k,to=200,by=1) 
  vec.TajimaD.max = f.TajimaD.max(n,vec.S,j,k)
  vec.TajimaD.min = f.TajimaD.min(n,vec.S,j,k)
  plot(vec.S,vec.TajimaD.max, xlab="No. of Segregating Sites", ylab="Tajima's D",col="orange",xlim=range(vec.S), ylim=range(c(vec.TajimaD.max, vec.TajimaD.min)),main=paste("k =",k,", n =",n),pch="o")
  lines(vec.S, vec.TajimaD.max, xlim=range(vec.S), ylim=range(vec.TajimaD.max), col="orange4",lwd=2,pch=16)
  points(vec.S,vec.TajimaD.min,col="pink",pch="o") 
  lines(vec.S, vec.TajimaD.min, xlim=range(vec.S), ylim=range(vec.TajimaD.min), col="pink4",lwd=2,pch=16)
  
}
# plot_tajima_D_S(n,S,j,k)

## Let's plot dependence of
## Tajima's D on n ##
# n.lower=3
# n.upper=50

#' Helper Function 
#' 
#' Comment 
#' 
#' @param n.lower comment
#' @param n.upper comment
#' @param S comment
#' @param j comment
#' @param k comment
#' @return the function
#' @export
#' 
plot_tajima_D_n <- function(n.lower, n.upper, S,j,k){
  vec.n = c(n.lower:n.upper)
  vec.TajimaD.max=c(1:length(vec.n))
  vec.TajimaD.min=c(1:length(vec.n))
  for (i in 1:length(vec.n)) {
    vec.TajimaD.max[i]=f.TajimaD.max(i+n.lower-1,S,j,k)
    vec.TajimaD.min[i]=f.TajimaD.min(i+n.lower-1,S,j,k) 
  }
  plot(vec.n,vec.TajimaD.max, xlab="No. of Sampled Individuals", ylab="Tajima's D",col="purple",xlim=range(vec.n), ylim=range(c(vec.TajimaD.max[!is.na(vec.TajimaD.max)], vec.TajimaD.min[!is.na(vec.TajimaD.max)])),main=paste("k =",k,", S =",S),pch="o")
  lines(vec.n, vec.TajimaD.max, xlim=range(vec.n), ylim=range(vec.TajimaD.max), col="darkslateblue",lwd=2,pch=16)

  points(vec.n,vec.TajimaD.min,col="green",pch="o") 
  lines(vec.n, vec.TajimaD.min, xlim=range(vec.n), ylim=range(vec.TajimaD.min), col="darkgreen",lwd=2,pch=16)
}

# plot_tajima_D_n(n.lower, n.upper, S,j,k)


#### Fu and Li's F ####

### Upper Bound ###
#' Helper Function 
#' 
#' Comment
#' 
#' @param n comment
#' @param S comment
#' @param j comment
#' @param k comment
#' @return the function
#' @export
#' 
f.FuLiF.max <- function(n,S,j,k) {
  (f.thetaPi.max(n,S,j,k) - k) / sqrt(u.F(n)*S + v.F(n)*S*S)
}

### Lower Bound ###
#' Helper Function 
#' 
#' Comment
#' 
#' @param n comment
#' @param S comment
#' @param j comment
#' @param k comment
#' @return the function
#' @export
#' 
f.FuLiF.min <- function(n,S,j,k) {
  (f.thetaPi.min(n,S,j,k) - k) / sqrt(u.F(n)*S + v.F(n)*S*S)
}

## Let's plot dependence of
## Fu and Li's F on frequency of singletons \xi_1 ##
#' Helper Function 
#' 
#' Comment
#' 
#' @param R comment
#' @param n comment
#' @param S comment
#' @param j comment
#' @param theta comment
#' @param scenario comment
#' @param plot_points comment
#' @param ... comment
#' @return the function
#' @export
#' 
plot_Fu_li_k <- function(R,n,S,j,theta=5,scenario ="neutral", plot_points = FALSE,...){
  vec.k = c(0:S)
  vec.FuLiF.max = f.FuLiF.max(n,S,j,vec.k)
  vec.FuLiF.min = f.FuLiF.min(n,S,j,vec.k)   
  plot(vec.k/S,vec.FuLiF.max, xlab="Frequency of Derived Singletons", ylab="Fu and Li's F",col="blue",
       type='l',lwd=2,ylim=range(c(vec.FuLiF.min, vec.FuLiF.max)),main=paste("n =",n,", S =",S),pch="o")
  #lines(vec.k/S, vec.FuLiF.max, xlim=range(vec.k/S), ylim=range(vec.FuLiF.max), col="blue",
        #lwd=2,pch=16)
  #points(vec.k/S,vec.FuLiF.min,col="red",pch="o") 
  lines(vec.k/S, vec.FuLiF.min, xlim=range(vec.k/S), ylim=range(vec.FuLiF.min), col="red",lwd=2,pch=16)
  if(plot_points == TRUE){
    runCommand(R,n,S,theta,f = 'FuLi',scenario)
    df = read.table("output.txt",header=FALSE)
    points(as.vector(df$V1[-1]),as.vector(df$V2[-1]),col="black",pch="o")
  }
}

## Let's plot dependence of
## Fu and Li's F on S ##
#' Helper Function 
#' 
#' Comment
#' 
#' @param n comment
#' @param j comment
#' @param k comment
#' @param ... comment
#' @return the function
#' @export
#' 
plot_Fu_li_S <- function(n,j,k,...){
  vec.S = seq(from=k,to=200,by=1) 
  vec.FuLiF.max = f.FuLiF.max(n,vec.S,j,k)
  vec.FuLiF.min = f.FuLiF.min(n,vec.S,j,k)
  plot(vec.S,vec.FuLiF.max, xlab="No. of Segregating Sites", ylab="Fu and Li's F",col="orange",xlim=range(vec.S), ylim=range(c(vec.FuLiF.min, vec.FuLiF.max)),main=paste("k =",k,", n =",n),pch="o")
  lines(vec.S, vec.FuLiF.max, xlim=range(vec.S), ylim=range(vec.FuLiF.max), col="orange4",lwd=2,pch=16)
  points(vec.S,vec.FuLiF.min,col="pink",pch="o") 
  lines(vec.S, vec.FuLiF.min, xlim=range(vec.S), ylim=range(vec.FuLiF.min), col="pink4",lwd=2,pch=16)
  
}

# plot_Fu_li_S(n,S,j,k)

## Let's plot dependence of
## Fu and Li's F on n ##

#' Helper Function 
#' 
#' Comment
#' 
#' @param n.lower comment
#' @param n.upper comment
#' @param S comment
#' @param j comment
#' @param k comment
#' @param ... comment
#' @return the function
#' @export
#' 
plot_Fu_li_n <- function(n.lower,n.upper,S,j,k,...){
  # n.lower=3
  # n.upper=100
  vec.n = c(n.lower:n.upper)
  vec.FuLiF.max=c(1:length(vec.n))
  vec.FuLiF.min=c(1:length(vec.n))
  for (i in 1:length(vec.n)) {
    vec.FuLiF.max[i]=f.FuLiF.max(i+n.lower-1,S,j,k)
    vec.FuLiF.min[i]=f.FuLiF.min(i+n.lower-1,S,j,k) 
  }
  plot(vec.n,vec.FuLiF.max, xlab="No. of Sampled Individuals", ylab="Fu and Li's F",col="purple",xlim=range(vec.n), ylim=range(c(vec.FuLiF.min, vec.FuLiF.max)),main=paste("k =",k,", S =",S),pch="o")
  lines(vec.n, vec.FuLiF.max, xlim=range(vec.n), ylim=range(vec.FuLiF.max), col="darkslateblue",lwd=2,pch=16)
  points(vec.n,vec.FuLiF.min,col="green",pch="o") 
  lines(vec.n, vec.FuLiF.min, xlim=range(vec.n), ylim=range(vec.FuLiF.min), col="darkgreen",lwd=2,pch=16)
}

# plot_Fu_li_n(n.lower,n.upper,S,j,k)


######### MEGA FUNCTIONS #######
#=============================================
#' A function to generate and plot summary statistics of parametric molecular genetic data,
#' together with their model-independent bounds.
#'
#' Return a plot specified by the input parameters.
#'
#' @param f a string denoting statistic {FuLi, TajPi, TajD}
#' @param R the number of data points to generate from simulation runs 
#' @param n the number of sampled individuals
#' @param S the number of segregating sites or polymorphic sites
#' @param j the argument of the component of the site frequency spectrum 
#' @param theta the mutation rate, need to be specified if plot_points=TRUE
#' @param scenario the parametric model {neutral, bottleneck, expansion}, if plot_points=TRUE
#' @param plot_points option of generating data from parametric model {TRUE, FALSE}
#' @param ... additional arguments to be passed (none in actuality)
#' @return plot of bounds according to parameter inputs, and optionally data points
#'
#' @keywords plot, simulation
#'
#' @export
#'
#' @examples
#' plot_statistic(f="FuLi",R=100,n=500,S=100,j=1,theta=5.0,scenario="neutral",plot_points=FALSE)
#' 

plot_statistic <- function(f = 'FuLi', R ,n, S,j, theta , scenario = 'neutral',
                           plot_points = FALSE, ...){
  if(f == 'FuLi'){
    plot_Fu_li_k(R,n,S,j,theta,scenario,plot_points)
  } else if (f == 'TajPi'){
    plot_tajima_Pi_k(R , n , S , j , theta,scenario,plot_points)
  } else if(f == 'TajD'){
    plot_tajima_D_k(R,n,S,j,theta,scenario,plot_points)
  }
}

# plot_statistic('TajPi', R = 100 ,n = 500,S = 100,j = 1, theta=5, scenario = 'neutral',plot_points = FALSE)
#=============================================

#' A function to plot model-independent bounds on summary statistics
#' as they depend on the number of individuals sampled.
#'
#' Return a plot specified by the input parameters.
#'
#' @param f a string denoting statistic {FuLi, TajPi, TajD}
#' @param n.lower the minimum number of sampled individuals, always set >=3
#' @param n.upper the maximum number of sampled individuals
#' @param S the number of segregating sites or polymorphic sites
#' @param j the argument of the component of the site frequency spectrum 
#' @param k the component of the site frequency spectrum
#' @param ... additional arguments to be passed (none in actuality)
#' @return plot of bounds according to parameter inputs
#'
#' @keywords plot, sampled individuals
#'
#' @export
#'
#' @examples
#' plot_n(f="FuLi",n.lower=3,n.upper=50,S=100,j=1,k=5)
#'

plot_n <- function(f= "FuLi",n.lower,n.upper, S,j,k, ...){
  if(f == 'FuLi'){
    plot_Fu_li_n(n.lower,n.upper,S,j,k)
  } else if (f == 'TajPi'){
    plot_tajima_Pi_n(n.lower,n.upper,S,j,k)
  } else if(f == 'TajD'){
    plot_tajima_D_n(n.lower,n.upper,S,j,k)
  }
}

#=============================================

#' A function to plot model-independent bounds on summary statistics
#' as they depend on the number of segregating or polymorphic sites.
#'
#' Return a plot specified by the input parameters.
#'
#' @param f a string denoting statistic {FuLi, TajPi, TajD}
#' @param n the number of sampled individuals
#' @param j the argument of the component of the site frequency spectrum 
#' @param k the component of the site frequency spectrum
#' @param ... additional arguments to be passed (none in actuality)
#' @return plot of bounds according to parameter inputs
#'
#' @keywords plot, sampled individuals
#'
#' @export
#'
#' @examples
#' plot_S(f="FuLi",n=50,j=1,k=5)
#'

plot_S <- function(f = 'FuLi',n,j,k,...){
  if(f == 'FuLi'){
    plot_Fu_li_S(n,j,k)
  } else if (f == 'TajPi'){
    plot_tajima_Pi_S(n,j,k)
  } else if(f == 'TajD'){
    plot_tajima_D_S(n,j,k)
  }
}

#=============================================

####### R code for WF Mutation Model #######
####### Version 2.0: (8/30/2017)     #######
####### By Alan Aw                   #######

### Setting Up ###
library(phyclust)
library(seqinr)

setwd("Downloads/msa")

### Helper Functions ###
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

## For Tajima's D ## (see Simonsen, Churchill, Aquadro, 1995)
v.T <- function(n) {
  (2*(n^2+n+3)/(9*n*(n-1))-(n+2)/(little.a(n)*n) + big.a(n)/(little.a(n))^2) / (little.a(n)^2+big.a(n))
}

u.T <- function(n) {
  (((n+1)/(3*(n-1))-1/little.a(n))/little.a(n)) - v.T(n) 
}

## For Fu and Li's F ##
v.F <- function(n) {
  ((2*n^3 + 110*n^2 -255*n + 153)/(9*n^2*(n-1))+(2*(n-1)*little.a(n))/(n^2)-(8*big.a(n))/n) / (little.a(n)^2+big.a(n))
}

u.F <- function(n) {
  ((4*n^2+19*n+3-12*(n+1)*little.a(n+1))/(3*n*(n-1)))/little.a(n) - v.F(n)
}

### Pre-simulation step ###
## HOW MANY REPS DO YOU WANT TO SIMULATE? ##
R = 1000
## (Optional) LET'S CREATE A DATA FRAME TO STORE SIMULATION OUTPUT! ## 
# df = data.frame(matrix(ncol=2,nrow=R,dimnames=list(c(),c("k","diffStat"))),stringsAsFactors=F)
k.vec = c()
diffStat.vec = c()


getCommand <- function(n, S, theta) {
  result = paste("./ms", n, "1 -t", theta, "-s", S, '-eN 0.2',theta, '> sfs.txt', sep = " ")
  result
}

runCommand <- function(R, command, S, n){
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
    #iOmega.vec[i] = (n-i)*i/choose(n,2) - 1/(little.a(n)) # Tajima's D (TOGGLE; uncomment to select) 
    iOmega.vec[i] = ifelse(i==1,1,0)*(2/n-1) + ifelse(i==1,0,1)*((n-i)*i/choose(n,2)) # Fu and Li's F (TOGGLE; uncomment to select)
  }
  # tajimaD=sum(iOmega.vec*sfs.vec)/sqrt(u.T(n)*S + v.T(n)*S*S) # Tajima's D (TOGGLE; uncomment to select)
  fuliF=sum(iOmega.vec*sfs.vec)/sqrt(u.F(n)*S + v.F(n)*S*S) # Fu and Li's F (TOGGLE; uncomment to select)
  ### This step adds simulated data to output list ###
  k.vec = c(k.vec,kOverS)
  # diffStat.vec=c(diffStat.vec,tajimaD) # Tajima's D (TOGGLE; uncomment to select)
  diffStat.vec=c(diffStat.vec,fuliF) # Fu and Li's F (TOGGLE; uncomment to select)

}

### This step produces the required data frame ###
## Specify simulation parameters ##
n = 500 # equals to first numeric argument in system formula 
theta = 5.0 # equals to "-t" in system formula
S = 100 # equals to "-s" in system formula

for (r in 1:R) {
  # system("./ms 500 1 -t 5.0 -s 100 > sfs.txt") # neutral (TOGGLE; uncomment to select)
  # system("./ms 500 1 -t 5.0 -s 100 -G 6.93 -eG 0.2 0.0 > sfs.txt") # population expansion (TOGGLE; uncomment to select)
  system("./ms 500 1 -t 5.0 -s 100 -eN 0.2 5.0 > sfs.txt") # population bottleneck (TOGGLE; uncomment to select)
  system("sed '1,7d' sfs.txt > sfs-modified.txt")
  #system("awk 'BEGIN{FS=" ";OFS=","}{$1=$1;print}' sfs-modified.txt > sfs.csv")
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
    #iOmega.vec[i] = (n-i)*i/choose(n,2) - 1/(little.a(n)) # Tajima's D (TOGGLE; uncomment to select) 
    iOmega.vec[i] = ifelse(i==1,1,0)*(2/n-1) + ifelse(i==1,0,1)*((n-i)*i/choose(n,2)) # Fu and Li's F (TOGGLE; uncomment to select)
  }
  # tajimaD=sum(iOmega.vec*sfs.vec)/sqrt(u.T(n)*S + v.T(n)*S*S) # Tajima's D (TOGGLE; uncomment to select)
  fuliF=sum(iOmega.vec*sfs.vec)/sqrt(u.F(n)*S + v.F(n)*S*S) # Fu and Li's F (TOGGLE; uncomment to select)
  ### This step adds simulated data to output list ###
  k.vec = c(k.vec,kOverS)
  # diffStat.vec=c(diffStat.vec,tajimaD) # Tajima's D (TOGGLE; uncomment to select)
  diffStat.vec=c(diffStat.vec,fuliF) # Fu and Li's F (TOGGLE; uncomment to select)
}

### This step produces the output .txt file of the simulation ###
simulation.output = data.frame(k.vec,diffStat.vec)
write.table(simulation.output, file = "fuliF-bottleneckN500S100.txt", row.names = FALSE)

################## END OF DATA SIMULATION! CONGRATULATIONS! ##################

## Example 2: Same as above, but with piping of output to sample_stats C program to generate summary statistics of dataset
system("./ms 20 1 -t 5.0 -s 1000 | ./sample_stats > samp.txt") # output file is .txt, can be changed to .ss, .csv, etc. 



# read table
  
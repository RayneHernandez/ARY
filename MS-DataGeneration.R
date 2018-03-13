n=500
S=100
R = 10
## (Optional) LET'S CREATE A DATA FRAME TO STORE SIMULATION OUTPUT! ## 
# df = data.frame(matrix(ncol=2,nrow=R,dimnames=list(c(),c("k","diffStat"))),stringsAsFactors=F)
k.vec = c()
diffStat.vec = c()

#system("cd Documents/GitProjects/ARY/")
for (r in 1:R) {
  # system("./ms 500 1 -t 5.0 -s 100 > sfs.txt") # neutral (TOGGLE; uncomment to select)
  # system("./ms 500 1 -t 5.0 -s 100 -G 6.93 -eG 0.2 0.0 > sfs.txt") # population expansion (TOGGLE; uncomment to select)
  working_dir = getwd()
  #result = paste('cd', working_dir)
  #result = paste(result,"&& ./ms 500 1 -t 5.0 -s 100 -eN 0.2 5.0 > sfs.txt")
  result = "./ms 500 1 -t 5.0 -s 100 -eN 0.2 5.0 > sfs.txt"
  system(result)
  #system(paste("cd",getwd() "&& ./ms 500 1 -t 5.0 -s 100 -eN 0.2 5.0 > sfs.txt",sep=" ")) # population bottleneck (TOGGLE; uncomment to select)
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
    iOmega.vec[i] = i*(n-i) / choose(n,2)
  }
  tajimaPi=sum(iOmega.vec*sfs.vec) # Tajima's Pi (TOGGLE; uncomment to select)
  ### This step adds simulated data to output list ###
  k.vec = c(k.vec,kOverS)
  # diffStat.vec=c(diffStat.vec,tajimaD) # Tajima's D (TOGGLE; uncomment to select)
  diffStat.vec=c(diffStat.vec,tajimaPi) # Fu and Li's F (TOGGLE; uncomment to select)
}

### This step produces the output .txt file of the simulation ###
simulation.output = data.frame(k.vec,diffStat.vec)
write.table(simulation.output, file = "tajimaPi-bottleneckN500S100.txt", row.names = FALSE)

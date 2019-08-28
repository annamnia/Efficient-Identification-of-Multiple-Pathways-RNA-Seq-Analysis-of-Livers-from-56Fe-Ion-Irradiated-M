##Written by Tianlong Chen, Sep, 2018
## unipartite network community detection, using modularity
## weighted or unweighted
#code to call the modularity algorithm, which was written in C.

# best power number picked by WGCNA, will change with networks 
pow <- 16   
#threshold for pearson correlation
thres <- 0.7  
#read in adjacency networks 
net <- read.table('adjacency_networks/adjacency_1e5.txt', header = T, row.names = 1)
net <- as.matrix(net)
#times to repeat running modularity algorithms
run <- 10
#size of networks
n <- dim(net)[1]

#Load community detection code into R
dyn.load("TC_ComDet_uni/R_uni.so") #computer specific

#Assign community detection code to a function in R
uni_cd <- function(size, net, iter, MODULARITY,ASSIGN){
  .C("uni_ComDet",as.integer(size),as.double(net),as.integer(iter),as.double(0.0),as.integer(1:size))}

#apply the threshold of pearson correlation
#set all edges with weight below the pearson correlation values to 0
for (i in 1:n) {
  net[i,i] <- 0
  for (j in 1:n) {
    if (net[i,j] < thres^pow) {
      net[i,j] <- 0
    }
  }
}

#removing disconnected genes
deg <- apply(net, 1, sum)
net <- net[deg != 0, deg != 0]
#network size
n <- dim(net)[1]

#run the modularity code on the data
modInput <- as.numeric(net)
modRes <- uni_cd(n, modInput, run)

#print network size
print(n)
#print modularity score
print(modRes[[4]])
#print number of clusters
print(max(modRes[[5]]))
#print median of cluster size
print(median(table(modRes[[5]])))
#print average cluster size
print(mean(table(modRes[[5]])))

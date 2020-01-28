## Written by Tianlong Chen, Sep, 2018
## unipartite network community detection, using modularity
## weighted or unweighted
## code to call the modularity algorithm, which was written in C.

library(svDialogs)

setwd(dlg_dir(default = getwd(), "Select Directory of Modularity Library", gui = .GUI)$res)
input_file<-dlg_open(default=getwd(), "Select Adjancect Matrix", multiple = FALSE, filters = dlg_filters["All", ],gui = .GUI)$res
output_file<-paste(tools::file_path_sans_ext(input_file),"modularity.csv",sep = "_")
output_file<-dlg_save(output_file, "Save Modularity Network", filters = dlg_filters["All", ],gui = .GUI)$res


# best power number picked by WGCNA, will change with networks 
pow<-as.numeric(dlg_input(message = "Enter Power Number from WCGNA", default = "16", gui = .GUI)$res)   
#threshold for pearson correlation
thres <-as.numeric(dlg_input(message = "Pearson Correlation Threshold", default = "0.7", gui = .GUI)$res)   
#read in adjacency networks 
net <- read.table(input_file[1], header = T, row.names = 1)
net <- as.matrix(net)
#times to repeat running modularity algorithms
run <-as.numeric(dlg_input(message = "Number of Times to Run Modularity Algorithm", default = "10", gui = .GUI)$res) 
#size of networks
n <- dim(net)[1]

#Load community detection code into R
if(.Platform$OS.typ=="windows")
{
  dyn.load("R_uni.dll") #windows
}
if(.Platform$OS.typ=="unix")
{
  dyn.load("R_uni.so") #unix
}

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
gene_assignments<-data.frame(as.vector(rownames(net)),modRes[5])
colnames(gene_assignments)<-c("Gene", "Cluster ID")
write.csv(gene_assignments,output_file)


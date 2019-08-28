
rm(list=ls())

# detach("package:convert", unload = TRUE)
# detach("package:Biobase", unload = TRUE)
# detach("package:edgeR", unload = TRUE)
# detach("package:DESeq2", unload = TRUE)
# detach("package:limma", unload = TRUE)
#2^0.59=1.51

source("https://bioconductor.org/biocLite.R")

biocLite(c("Biobase"))
biocLite("convert")
biocLite("edgeR")
biocLite("GenomicRanges")
a
biocLite("limma")
install.packages("xlsx")
install.packages("rJava")
install.packages("glue")


library(dynamicTreeCut)
library(Biobase)
library(convert)
library(IRanges)
library(edgeR)
library(GenomicRanges)
library(DESeq2)
library(limma)
library(rJava)
library(xlsx)
library(data.table)
library(WGCNA)
library(dplyr)
library(readr)
library(tidyr)


#cutreeDynamic
#deepSplit(0,1,2,3,4) controls how finely clusters will be split
#pamStage (FALSE o TRUE) truns PAM stage off/on 



sessionInfo()

setwd("C:/Users/amnia/Dropbox/IPA_Paper/WGCNA")

exprsFile <- "exprsData.txt"

exprs <- as.matrix(read.table(exprsFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE))
pDataFile<- "pData2.txt"
pData <-read.table(pDataFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE)
all(rownames(pData)==colnames(exprs))
# sapply(pData,class)
# pData[c(15,20),c("Strain","Treatment")]
metadata<- data.frame(labelDescription= c("Mouse Strain","Radiated/Non_Iradiated", "Age of Sacrifice"),
                      row.names = c("Strain","Treatment","Time"))

phenoData<- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
head(pData(phenoData))
phenoData[c("B1","B2"),"Treatment"]
exampleSet<- ExpressionSet(assayData=exprs,phenoData=phenoData)

# Plot the entire dataset
# cpm.mat <- log(cpm(exprs(exampleSet)))
# mean.vec <- apply(cpm.mat, 1, mean)
# sdvec <- apply(cpm.mat, 1, sd)
# plot(mean.vec, sdvec, pch=".", main="3 Replicates", ylab="sd", xlab="Average logCPM")
# dev.copy(png,'Entire_Data.png')
# dev.off()



#Create edgeR object
#edgeR uses TMM normalization method
# Create the contrast matrix
#if you don't have zero you're just comparing one against the other, C57 or not C57, it does it alphabetically
#one group vs the other group
## edgeR ##
# Design matrix
#The nested interaction model makes it easy to nd genes that respond to the treatment at any
# time, in a single test. Continuing the above example,
# nds genes that respond to the treatment at either 1 hour or 2 hours versus the 0 hour baseline.
# This is analogous to an ANOVA F-test for a normal linear model.
# design <- model.matrix(~0+Group)
# colnames(design) <- levels(Group)

pData(exampleSet)$Treatment<-as.factor(pData(exampleSet)$Treatment)
pData(exampleSet)$Time<- as.factor(pData(exampleSet)$Time)
pData(exampleSet)$Strain<- as.factor(pData(exampleSet)$Strain)

pData(exampleSet)$Treatment <- relevel(pData(exampleSet)$Treatment, ref="Control")
pData(exampleSet)$Time<-relevel(pData(exampleSet)$Time, ref="1mo")
pData(exampleSet)$Strain<-relevel(pData(exampleSet)$Strain, ref="C57")

#set threshold
p.threshold <- 0.05
Group<-factor(paste(pData(exampleSet)$Strain,pData(exampleSet)$Treatment,pData(exampleSet)$Time,sep="."))
cbind(pData(exampleSet),Group=Group)

design <- model.matrix(~Treatment, data=pData(exampleSet))
length(colnames(design))
View(colnames(design))

# design <- model.matrix(~0+Group, data=pData(exampleSet))
# length(colnames(design))
# View(colnames(design))

#7.9
dge <- DGEList(counts=exprs(exampleSet), group=Group)
#filter out lowly expressed genes 
keep <- rowSums(cpm(dge)>1) >= 2
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge.table<-dge$counts
# write.xlsx(dge.table, "C:/Users/amnia/Dropbox/IPA_Paper/WGCNA/edgeR_Filtered.xlsx",
#            col.names=TRUE, row.names=TRUE, append=FALSE)

write.csv(dge.table,file="edgeR_Filtered.csv")


#DESeq
exprsFile <- "edgeR_Filtered.txt"
exprs <- as.matrix(read.table(exprsFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE))
pDataFile<- "pData2.txt"
pData <-read.table(pDataFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE)
all(rownames(pData)==colnames(exprs))
# sapply(pData,class)
# pData[c(15,20),c("Strain","Treatment")]
metadata<- data.frame(labelDescription= c("Mouse Strain","Radiated/Non_Iradiated", "Age of Sacrifice"),
                      row.names = c("Strain","Treatment","Time"))

phenoData<- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
# head(pData(phenoData))
# phenoData[c("B1","B2"),"Treatment"]
exampleSet<- ExpressionSet(assayData=exprs,phenoData=phenoData)

pData(exampleSet)$Treatment<-as.factor(pData(exampleSet)$Treatment)
pData(exampleSet)$Time<- as.factor(pData(exampleSet)$Time)
pData(exampleSet)$Strain<- as.factor(pData(exampleSet)$Strain)

pData(exampleSet)$Treatment <- relevel(pData(exampleSet)$Treatment, ref="Control")
pData(exampleSet)$Time<-relevel(pData(exampleSet)$Time, ref="1mo")
pData(exampleSet)$Strain<-relevel(pData(exampleSet)$Strain, ref="C57")
#set threshold
p.threshold <- 0.05
Group<-factor(paste(pData(exampleSet)$Strain,pData(exampleSet)$Treatment,pData(exampleSet)$Time,sep="."))
#cbind(pData(exampleSet),Group=Group)
design <- model.matrix(~Treatment, data=pData(exampleSet))
length(colnames(design))
View(colnames(design))

#create DESeq2 dataset
dds<-DESeqDataSetFromMatrix(countData=exprs(exampleSet), colData=pData(exampleSet), design=~Treatment)
dds <- DESeq(dds)
#varianceStabilzingTransformtion
varianceStabilizingTransformation(dds, blind = TRUE,
                                  fitType = "parametric")

deseq_normalized<-getVarianceStabilizedData(dds)
# write.xlsx(deseq_normalized, "C:/Users/amnia/Dropbox/IPA_Paper/WGCNA/deseq_normalized.xlsx",
#            col.names=TRUE, row.names=TRUE, append=FALSE)


#WGCNA
write.csv(deseq_normalized,file="deseq_normalized.csv")


exprData = read.csv("deseq_normalized.csv",stringsAsFactors = FALSE)
dim(exprData)
names(exprData)

datExpr0 = as.data.frame(t(exprData[,-c(1)]))
dim(datExpr0)
names(datExpr0) = exprData$Gene_Ensemble_ID
rownames(datExpr0) = names(exprData)[-c(1)]
rownames(datExpr0)

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

dev.off()
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(6,5)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


traitData = read.delim("pData3.txt");
dim(traitData)
names(traitData)

# Form a data frame analogous to expression data that will hold the clinical traits.
Samples = rownames(datExpr0)
traitRows = match(Samples, traitData$Sample)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]
collectGarbage();


dev.off()
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

dev.off()
#Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 17, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



dev.off()
#PNAS paper power
powers1=c(seq(1,10,by=1),seq(73,83,by=5))
RpowerTable = pickSoftThreshold(datExpr0, powerVector= powers1,RsquaredCut = 0.90)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(RpowerTable$fitIndices[,1], -sign(RpowerTable$fitIndices[,3])*RpowerTable$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(RpowerTable$fitIndices[,1], -sign(RpowerTable$fitIndices[,3])*RpowerTable$fitIndices[,2],
     labels=powers1,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(RpowerTable$fitIndices[,1], RpowerTable$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(RpowerTable$fitIndices[,1], RpowerTable$fitIndices[,5], labels=powers1, cex=cex1,col="red")





softPower = 73
#Defaults to Pearson correlation.
adjacency = adjacency(datExpr0, power = softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

#added on 09/06/2018
sample <- dissTOM[c(1:10),c(1:10)]
library(igraph)
g <- graph_from_adjacency_matrix(data.matrix(sample), mode = 'undirected', weighted = T)
plot(g, layout = layout_with_kk)


write.table(TOM,file="TOM.txt")
write.csv(TOM,file="TOM.csv")

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "Fe_C57_124Mo_networkConstruction-stepByStep.RData")

# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

red.module=names(datExpr0)[moduleColors=="red"]
write.xlsx(red.module, "C:/Users/amnia/Dropbox/IPA_Paper/WGCNA/red.module.xlsx",
           col.names=TRUE, row.names=TRUE, append=FALSE)


turquoise.module=names(datExpr0)[moduleColors=="turquoise"]
write.xlsx(turquoise.module, "C:/Users/amnia/Dropbox/IPA_Paper/WGCNA/turquoise.module.xlsx",
           col.names=TRUE, row.names=TRUE, append=FALSE)


black.module=names(datExpr0)[moduleColors=="black"]
write.xlsx(black.module, "C:/Users/amnia/Dropbox/IPA_Paper/WGCNA/black.module.xlsx",
           col.names=TRUE, row.names=TRUE, append=FALSE)

yellow.module=names(datExpr0)[moduleColors=="yellow"]
write.xlsx(yellow.module, "C:/Users/amnia/Dropbox/IPA_Paper/WGCNA/yellow.module.xlsx",
           col.names=TRUE, row.names=TRUE, append=FALSE)

# Define variable weight containing the weight column of datTrait
Radiation_Status = as.data.frame(datTraits$Treatment);
names(Radiation_Status) = "Radiation_Status"

Time = as.data.frame(datTraits$Time);
names(Time) = "Time"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, Radiation_Status, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Radiation_Status), sep="");
names(GSPvalue) = paste("p.GS.", names(Radiation_Status), sep="");



# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, Time, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Time), sep="");
names(GSPvalue) = paste("p.GS.", names(Time), sep="");



dev.off()
module = "yellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


dev.off()
module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Time",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, Radiation_Status, use = "p")));


multExpr=list(datExpr0,(datExpr0)[moduleColors=="red"])
multcolor=list()
mp=modulePreservation()

# codetime<-Sys.Date()
# writeLines(capture.output(sessionInfo()),"sessionInfo.txt")
save.image(file ="wgcna_Anna_C57_Fe_alltime.RData")


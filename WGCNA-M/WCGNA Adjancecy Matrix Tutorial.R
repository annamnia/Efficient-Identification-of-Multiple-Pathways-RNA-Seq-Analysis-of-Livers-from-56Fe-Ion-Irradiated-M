rm(list=ls())
source("http://bioconductor.org/biocLite.R")

library(gdata)
library(dplyr)
library(Biobase)
library(convert)
library(IRanges)
library(edgeR)
library(GenomicRanges)
library(DESeq2)
library(limma)
library(rJava)
library(xlsx)
library(biomaRt)
library(enrichR)
library(devtools)
library(WGCNA)

#C57_Fe_1mo
sessionInfo()
path<-setwd("/Users/WGCNA_M/1_BMC_Bioinformatics_Methodology/Scripts/Modularity Example/") 
#change the working directory path to the location where all of the example files are located. 
p.threshold <- 0.05

exprsFile <- "raw_exprsData_Control_C57_1mo.txt"
exprs <- as.matrix(read.table(exprsFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE))
pDataFile<- "pData_Fe_Control_C57_1mo.txt"
pData <-read.table(pDataFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE)
all(rownames(pData)==colnames(exprs))
metadata<- data.frame(labelDescription= c("Mouse Strain","Radiated/Non_Iradiated", "Age of Sacrifice"),
                      row.names = c("Strain","Treatment","Time"))

phenoData<- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
exampleSet<- ExpressionSet(assayData=exprs,phenoData=phenoData)


pData(exampleSet)$Treatment<-as.factor(pData(exampleSet)$Treatment)
pData(exampleSet)$Time<- as.factor(pData(exampleSet)$Time)
pData(exampleSet)$Strain<- as.factor(pData(exampleSet)$Strain)

pData(exampleSet)$Treatment <- relevel(pData(exampleSet)$Treatment, ref="Control")

Group<-factor(paste(pData(exampleSet)$Strain,pData(exampleSet)$Treatment,pData(exampleSet)$Time,sep="."))
cbind(pData(exampleSet),Group=Group)
design <- model.matrix(~Treatment, data=pData(exampleSet))
length(colnames(design))

dge <- DGEList(counts=exprs(exampleSet), group=Group)
# Normalize by total count
dge <- calcNormFactors(dge)
# Estimate dispersion parameter for GLM
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design, method="power")
dge<- estimateGLMTagwiseDisp(dge,design)

# Model fitting
#quasi-likelihood negative binomial generalized log-linear model to count data
fit.edgeR <- glmQLFit(dge, design)
names(fit.edgeR)

lrt.edgeR <- glmLRT(fit.edgeR, coef =2)
edgeR_results <- lrt.edgeR$table
sig.edgeR <- decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value = p.threshold)
genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR != 0)]
topgenes<-topTags(lrt.edgeR, n=length(genes.edgeR), adjust.method="BH", sort.by="PValue")
keep<-topgenes$table$FDR <= p.threshold  & abs(topgenes$table$logFC) >= 0.59
Total_DGE<-dim(topgenes[keep,]$table)[[1]]
upreg<-sum(topgenes$table$logFC>= 0.59)
downreg<-sum(topgenes$table$logFC<=-0.59)
Summary<-c("Total DEG",Total_DGE,"Upregulated",upreg,"Downregulated",downreg)
write.table(Summary,paste(path,"/C57_1mo_TreatmentFe_Summary.txt",sep=""),
            col.names=FALSE, row.names=FALSE, append=FALSE)
write.xlsx(topgenes[keep,]$table, paste(path,"/edgeR_DGE_C57_1mo_TreatmentFe.xlsx",sep=""),
           col.names=TRUE, row.names=TRUE, append=FALSE)


#C57_Fe_2mo
sessionInfo()
exprsFile <- "raw_exprsData_Control_C57_2mo.txt"
exprs <- as.matrix(read.table(exprsFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE))
pDataFile<- "pData_Fe_Control_C57_2mo.txt"
pData <-read.table(pDataFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE)
all(rownames(pData)==colnames(exprs))

metadata<- data.frame(labelDescription= c("Mouse Strain","Radiated/Non_Iradiated", "Age of Sacrifice"),
                      row.names = c("Strain","Treatment","Time"))

phenoData<- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)

exampleSet<- ExpressionSet(assayData=exprs,phenoData=phenoData)

pData(exampleSet)$Treatment<-as.factor(pData(exampleSet)$Treatment)
pData(exampleSet)$Time<- as.factor(pData(exampleSet)$Time)
pData(exampleSet)$Strain<- as.factor(pData(exampleSet)$Strain)

pData(exampleSet)$Treatment <- relevel(pData(exampleSet)$Treatment, ref="Control")


Group<-factor(paste(pData(exampleSet)$Strain,pData(exampleSet)$Treatment,pData(exampleSet)$Time,sep="."))
cbind(pData(exampleSet),Group=Group)
design <- model.matrix(~Treatment, data=pData(exampleSet))
length(colnames(design))

dge <- DGEList(counts=exprs(exampleSet), group=Group)
# Normalize by total count
dge <- calcNormFactors(dge)
# Estimate dispersion parameter for GLM
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design, method="power")
dge<- estimateGLMTagwiseDisp(dge,design)

# Model fitting
#quasi-likelihood negative binomial generalized log-linear model to count data
fit.edgeR <- glmQLFit(dge, design)
names(fit.edgeR)

#Coef=2, TreatmentFe:2mo
lrt.edgeR <- glmLRT(fit.edgeR, coef =2)
edgeR_results <- lrt.edgeR$table
sig.edgeR <- decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value = p.threshold)
genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR != 0)]
topgenes<-topTags(lrt.edgeR, n=length(genes.edgeR), adjust.method="BH", sort.by="PValue")
keep<-topgenes$table$FDR <= p.threshold & abs(topgenes$table$logFC) >= 0.59
Total_DGE<-dim(topgenes[keep,]$table)[[1]]
upreg<-sum(topgenes$table$logFC>= 0.59)
downreg<-sum(topgenes$table$logFC<=-0.59)
Summary<-c("Total DEG",Total_DGE,"Upregulated",upreg,"Downregulated",downreg)
write.table(Summary,paste(path,"/C57_2mo_TreatmentFe_Summary.txt",sep=""),
            col.names=FALSE, row.names=FALSE, append=FALSE)
write.xlsx(topgenes[keep,]$table, paste(path,"/edgeR_DGE_C57_2mo_TreatmentFe.xlsx",sep=""),
           col.names=TRUE, row.names=TRUE, append=FALSE)

#C57_Fe_4mo
sessionInfo()
exprsFile <- "raw_exprsData_Control_C57_4mo.txt"
exprs <- as.matrix(read.table(exprsFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE))
pDataFile<- "pData_Fe_Control_C57_4mo.txt"
pData <-read.table(pDataFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE)
all(rownames(pData)==colnames(exprs))
# sapply(pData,class)
# pData[c(15,20),c("Strain","Treatment")]
metadata<- data.frame(labelDescription= c("Mouse Strain","Radiated/Non_Iradiated", "Age of Sacrifice"),
                      row.names = c("Strain","Treatment","Time"))

phenoData<- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
exampleSet<- ExpressionSet(assayData=exprs,phenoData=phenoData)

pData(exampleSet)$Treatment<-as.factor(pData(exampleSet)$Treatment)
pData(exampleSet)$Time<- as.factor(pData(exampleSet)$Time)
pData(exampleSet)$Strain<- as.factor(pData(exampleSet)$Strain)

pData(exampleSet)$Treatment <- relevel(pData(exampleSet)$Treatment, ref="Control")


#set threshold
p.threshold <- 0.05
Group<-factor(paste(pData(exampleSet)$Strain,pData(exampleSet)$Treatment,pData(exampleSet)$Time,sep="."))
cbind(pData(exampleSet),Group=Group)
design <- model.matrix(~Treatment, data=pData(exampleSet))
length(colnames(design))

dge <- DGEList(counts=exprs(exampleSet), group=Group)
# Normalize by total count
dge <- calcNormFactors(dge)
# Estimate dispersion parameter for GLM
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design, method="power")
dge<- estimateGLMTagwiseDisp(dge,design)

# Model fitting
#quasi-likelihood negative binomial generalized log-linear model to count data
fit.edgeR <- glmQLFit(dge, design)
names(fit.edgeR)

#Coef=2, TreatmentFe:4mo
lrt.edgeR <- glmLRT(fit.edgeR, coef =2)
edgeR_results <- lrt.edgeR$table
sig.edgeR <- decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value = p.threshold)
genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR != 0)]
topgenes<-topTags(lrt.edgeR, n=length(genes.edgeR), adjust.method="BH", sort.by="PValue")
keep<-topgenes$table$FDR <= p.threshold & abs(topgenes$table$logFC) >= 0.59
Total_DGE<-dim(topgenes[keep,]$table)[[1]]
upreg<-sum(topgenes$table$logFC>= 0.59)
downreg<-sum(topgenes$table$logFC<=-0.59)
Summary<-c("Total DEG",Total_DGE,"Upregulated",upreg,"Downregulated",downreg)
write.table(Summary,paste(path,"/C57_4mo_TreatmentFe_Summary.txt",sep=""),
            col.names=FALSE, row.names=FALSE, append=FALSE)
write.xlsx(topgenes[keep,]$table, paste(path,"/edgeR_DGE_C57_4mo_TreatmentFe.xlsx",sep=""),
           col.names=TRUE, row.names=TRUE, append=FALSE)



#C57_Fe_9mo
sessionInfo()
exprsFile <- "raw_exprsData_Control_C57_9mo.txt"
exprs <- as.matrix(read.table(exprsFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE))
pDataFile<- "pData_Fe_Control_C57_9mo.txt"
pData <-read.table(pDataFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE)
all(rownames(pData)==colnames(exprs))
# sapply(pData,class)
# pData[c(15,20),c("Strain","Treatment")]
metadata<- data.frame(labelDescription= c("Mouse Strain","Radiated/Non_Iradiated", "Age of Sacrifice"),
                      row.names = c("Strain","Treatment","Time"))

phenoData<- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
exampleSet<- ExpressionSet(assayData=exprs,phenoData=phenoData)

pData(exampleSet)$Treatment<-as.factor(pData(exampleSet)$Treatment)
pData(exampleSet)$Time<- as.factor(pData(exampleSet)$Time)
pData(exampleSet)$Strain<- as.factor(pData(exampleSet)$Strain)

pData(exampleSet)$Treatment <- relevel(pData(exampleSet)$Treatment, ref="Control")

Group<-factor(paste(pData(exampleSet)$Strain,pData(exampleSet)$Treatment,pData(exampleSet)$Time,sep="."))
cbind(pData(exampleSet),Group=Group)
design <- model.matrix(~Treatment, data=pData(exampleSet))
length(colnames(design))

dge <- DGEList(counts=exprs(exampleSet), group=Group)
# Normalize by total count
dge <- calcNormFactors(dge)
# Estimate dispersion parameter for GLM
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design, method="power")
dge<- estimateGLMTagwiseDisp(dge,design)

# Model fitting
#quasi-likelihood negative binomial generalized log-linear model to count data
fit.edgeR <- glmQLFit(dge, design)
names(fit.edgeR)

#Coef=2, TreatmentFe:9mo
lrt.edgeR <- glmLRT(fit.edgeR, coef =2)
edgeR_results <- lrt.edgeR$table
sig.edgeR <- decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value = p.threshold)
genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR != 0)]
topgenes<-topTags(lrt.edgeR, n=length(genes.edgeR), adjust.method="BH", sort.by="PValue")
keep<-topgenes$table$FDR <= p.threshold & abs(topgenes$table$logFC) >= 0.59
Total_DGE<-dim(topgenes[keep,]$table)[[1]]
upreg<-sum(topgenes$table$logFC>= 0.59)
downreg<-sum(topgenes$table$logFC<=-0.59)
Summary<-c("Total DEG",Total_DGE,"Upregulated",upreg,"Downregulated",downreg)
write.table(Summary,paste(path,"/C57_9mo_TreatmentFe_Summary.txt",sep=""),
            col.names=FALSE, row.names=FALSE, append=FALSE)
write.xlsx(topgenes[keep,]$table, paste(path,"/edgeR_DGE_C57_9mo_TreatmentFe.xlsx",sep=""),
           col.names=TRUE, row.names=TRUE, append=FALSE)


#C57_Fe_12mo
sessionInfo()
exprsFile <- "raw_exprsData_Control_C57_12mo.txt"
exprs <- as.matrix(read.table(exprsFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE))
pDataFile<- "pData_Fe_Control_C57_12mo.txt"
pData <-read.table(pDataFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE)
all(rownames(pData)==colnames(exprs))
metadata<- data.frame(labelDescription= c("Mouse Strain","Radiated/Non_Iradiated", "Age of Sacrifice"),
                      row.names = c("Strain","Treatment","Time"))

phenoData<- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
exampleSet<- ExpressionSet(assayData=exprs,phenoData=phenoData)

pData(exampleSet)$Treatment<-as.factor(pData(exampleSet)$Treatment)
pData(exampleSet)$Time<- as.factor(pData(exampleSet)$Time)
pData(exampleSet)$Strain<- as.factor(pData(exampleSet)$Strain)

pData(exampleSet)$Treatment <- relevel(pData(exampleSet)$Treatment, ref="Control")

Group<-factor(paste(pData(exampleSet)$Strain,pData(exampleSet)$Treatment,pData(exampleSet)$Time,sep="."))
cbind(pData(exampleSet),Group=Group)
design <- model.matrix(~Treatment, data=pData(exampleSet))
length(colnames(design))

dge <- DGEList(counts=exprs(exampleSet), group=Group)
# Normalize by total count
dge <- calcNormFactors(dge)
# Estimate dispersion parameter for GLM
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design, method="power")
dge<- estimateGLMTagwiseDisp(dge,design)

# Model fitting
#quasi-likelihood negative binomial generalized log-linear model to count data
fit.edgeR <- glmQLFit(dge, design)
names(fit.edgeR)

#Coef=2, TreatmentFe:12
lrt.edgeR <- glmLRT(fit.edgeR, coef =2)
edgeR_results <- lrt.edgeR$table
sig.edgeR <- decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value = p.threshold)
genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR != 0)]
topgenes<-topTags(lrt.edgeR, n=length(genes.edgeR), adjust.method="BH", sort.by="PValue")
keep<-topgenes$table$FDR <= p.threshold & abs(topgenes$table$logFC) >= 0.59
Total_DGE<-dim(topgenes[keep,]$table)[[1]]
upreg<-sum(topgenes$table$logFC>= 0.59)
downreg<-sum(topgenes$table$logFC<=-0.59)
Summary<-c("Total DEG",Total_DGE,"Upregulated",upreg,"Downregulated",downreg)
write.table(Summary,paste(path,"/C57_12mo_TreatmentFe_Summary.txt",sep=""),
            col.names=FALSE, row.names=FALSE, append=FALSE)
write.xlsx(topgenes[keep,]$table, paste(path,"/edgeR_DGE_C57_12mo_TreatmentFe.xlsx",sep=""),
           col.names=TRUE, row.names=TRUE, append=FALSE)

#now take the files from DGE and normalize the raw files with DESeq

C57_Fe_1mo <- read.xlsx("edgeR_DGE_C57_1mo_TreatmentFe.xlsx", sheetName = "Sheet1")
C57_Fe_2mo <- read.xlsx("edgeR_DGE_C57_2mo_TreatmentFe.xlsx", sheetName = "Sheet1")
C57_Fe_4mo <- read.xlsx("edgeR_DGE_C57_4mo_TreatmentFe.xlsx",sheetName = "Sheet1")
C57_Fe_9mo <- read.xlsx("edgeR_DGE_C57_9mo_TreatmentFe.xlsx",sheetName = "Sheet1")
C57_Fe_12mo <- read.xlsx("edgeR_DGE_C57_12mo_TreatmentFe.xlsx",sheetName = "Sheet1")

C57_Fe_alltime <- rbind(C57_Fe_1mo,C57_Fe_2mo,C57_Fe_4mo,C57_Fe_9mo,C57_Fe_12mo)
C57_Fe_alltime_tb<- tbl_df(C57_Fe_alltime)
View(C57_Fe_alltime_tb)
names(C57_Fe_alltime_tb)
colnames(C57_Fe_alltime_tb) <- c("Ensemble_Gene_ID","logFC","logCPM","LR","PValue","FDR")
C57_Fe_alltime_tb<-arrange(C57_Fe_alltime_tb,FDR)


z<-C57_Fe_alltime_tb[!duplicated(C57_Fe_alltime_tb$Ensemble_Gene_ID),]
z2<-C57_Fe_alltime_tb[duplicated(C57_Fe_alltime_tb$Ensemble_Gene_ID),]

write.csv(z,file="not_duplicated_.05_FDR.csv")
write.csv(z2,file="duplicated.csv")
write.csv(C57_Fe_alltime_tb,file="all.csv")


#head(C57_Fe_alltime_tb)

FDR_threshold=0.00001
not_duplicated_FDR<-z%>% filter(FDR<=FDR_threshold)
duplicated_FDR<-z2%>% filter(FDR<=FDR_threshold)
C57_Fe_alltime_filtered<-C57_Fe_alltime_tb%>% filter(FDR<=FDR_threshold)
write.csv(not_duplicated_FDR,file="not_duplicated_FDR_threshold.csv")
write.csv(duplicated_FDR,file="duplicated_FDR_threshold.csv")
write.csv(C57_Fe_alltime_filtered,file="all_FDR_threshold.csv")


#Make a DESeq normalized file from raw data files
exprsFile <- "raw_exprsData_C57.txt"
exprs <- as.matrix(read.table(exprsFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE))

pDataFile<- "pData_Fe_C57.txt"
pData <-read.table(pDataFile, row.names = 1, header=TRUE, sep="\t", as.is=TRUE)
all(rownames(pData)==colnames(exprs))
# # sapply(pData,class)
# # pData[c(15,20),c("Strain","Treatment")]
metadata<- data.frame(labelDescription= c("Mouse Strain","Radiated/Non_Iradiated", "Age of Sacrifice"),
                      row.names = c("Strain","Treatment","Time"))
phenoData<- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)

exampleSet<- ExpressionSet(assayData=exprs,phenoData=phenoData)
pData(exampleSet)$Treatment<-as.factor(pData(exampleSet)$Treatment)
pData(exampleSet)$Time<- as.factor(pData(exampleSet)$Time)
pData(exampleSet)$Strain<- as.factor(pData(exampleSet)$Strain)

pData(exampleSet)$Treatment <- relevel(pData(exampleSet)$Treatment, ref="Control")
pData(exampleSet)$Time<-relevel(pData(exampleSet)$Time, ref="1mo")
pData(exampleSet)$Strain<-relevel(pData(exampleSet)$Strain, ref="C57")



design <- model.matrix(~Treatment, data=pData(exampleSet))

# #create DESeq2 dataset
dds<-DESeqDataSetFromMatrix(countData=exprs(exampleSet), colData=pData(exampleSet), design=~Treatment)
dds <- DESeq(dds)
# #varianceStabilzingTransformtion
varianceStabilizingTransformation(dds, blind = TRUE,
                                  fitType = "parametric")

deseq_normalized<-getVarianceStabilizedData(dds)
write.csv(deseq_normalized,file="deseq_normalized.csv")


#retrieve the FDR<10^-5 genes nonduplicated from the deseq normalized file and read into an exprs file for WGCNA
not_duplicated<- read.csv("not_duplicated_FDR_threshold.csv", stringsAsFactors = FALSE)
deseq_normalized<- read.csv("deseq_normalized.csv", stringsAsFactors = FALSE)

#match gives relative position of the first one in the second
m<-match(not_duplicated$Ensemble_Gene_ID,deseq_normalized$X)
exprsData<-deseq_normalized[m,]
write.csv(exprsData,file ="exprsData.csv")

#WGCNA Adjacency Calculation
exprData = read.csv("exprsData.csv",stringsAsFactors = FALSE)
dim(exprData)
names(exprData)


datExpr0 = as.data.frame(t(exprData[,-c(1,2)]))
dim(datExpr0)
names(datExpr0) = exprData$X
rownames(datExpr0) = names(exprData)[-c(1,2)]
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
collectGarbage()


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
#tutorial
#Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
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

softPower = 16
#Defaults to Pearson correlation.
adjacency = adjacency(datExpr0, power = softPower)
#write adj into a csv file and go to mod or dynamic tree cut
write.table(adjacency,"adjancency_matrix.txt",sep='\t')

#two paths: continue with dynamic tree cut or go to modularity


codetime<-Sys.Date()
writeLines(capture.output(sessionInfo()),"sessionInfo.txt")
save.image(file =".RData")
#the document was processed on :





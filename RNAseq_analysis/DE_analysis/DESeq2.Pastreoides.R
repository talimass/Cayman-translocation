# necessary libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(ggforce)
library(VennDiagram)

# setting working directory 
setwd("/home/gospozha/haifa/cayman/rna/mapping/github")

# listing files with gene counts for each sample
dir = "/home/gospozha/haifa/cayman/rna/mapping/"
files = list.files(paste0(dir, "count.tables"), "*ReadsPerGene.out.tab", full.names = T)
countData = data.frame(fread(files[1]))[c(1,2)]

# looping and reading the 2nd column from the remaining files
for(i in 2:length(files)) {
  countData = cbind(countData, data.frame(fread(files[i]))[2])
}

# skipping the first 4 lines, since count data starts on the 5th line
countData = countData[c(5:nrow(countData)),]

# renaming columns as sample names
colnames(countData) = c("GeneID", gsub(paste0(dir,"count.tables/"), "", files))
colnames(countData) = gsub("_ReadsPerGene.out.tab", "", colnames(countData))
rownames(countData) = countData$GeneID
countData = countData[,c(2:ncol(countData))]
names <- colnames(countData)

# writing count matrix to a file
write.csv(countData, file="CountMatrix.csv")  

#### shallow ####

# reading count matrix from a file
countData  <- read.csv2('CountMatrix.csv', header=TRUE, row.names=1, sep=',')
colnames(countData) <- names
# reading metadata file
MetaData <- read.csv2('Metadata.origin.S.csv', header=TRUE, sep=',')
# removing the samples without any info in MetaData from the count table
countData <- countData[ , MetaData$id] 
# changing columns to factors
MetaData$condition <- factor(MetaData$condition, levels = c("S","SS","SD"))
MetaData$origin <- as.factor(MetaData$origin)
# creating a DESeq2 object 
# since site is already nested within the origin factor, we don't need to
# additionally account for it - origin already absorbs the site effect 
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = MetaData,
                              design = ~ origin  + condition)

# pre filtering - removing rows with low counts
smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# running a model
dds <- DESeq(dds)

# estimating size factors to determine if we can use vst
# to transform our data. Size factors should be less than 4 to use vst
SF <- estimateSizeFactors(dds) 
print(sizeFactors(SF))
# variance-stabilizing transformation
vst <- varianceStabilizingTransformation(dds)

# since vst does not remove variation that can be associated with covariates, 
# we manually remove the effect of covariates to be able to visualize it on PCA
mat <- assay(vst)
mm <- model.matrix(~condition, colData(vst))
mat <- limma::removeBatchEffect(mat, batch=vst$origin, design=mm)
assay(vst) <- mat
# PCA plot
pcaData <- plotPCA(vst, intgroup=c("condition", "site"), ntop = 1000, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("PCA.shallow.pdf",width=7)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=site)) +
  geom_point(size=3) +
  ggtitle("PCA plot for translocation from shallow reef") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()
dev.off()

# PCA with ellipses
pdf("PCA.shallow.ellips.pdf", width=7)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  geom_mark_ellipse(aes(fill=condition), alpha=0.3) +
  ggtitle("PCA plot for translocation from shallow reef") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()
dev.off()

# sample distances
sampleDists <- dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$condition, vst$site, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("Dist.shallow.pdf",width=7)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

# LRT test to test the effect of depth in general
# removing "depth" to test its effect
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ origin)  
res_lrt <- results(dds_lrt)
summary(res_lrt)
# 8.2 up, 2.4 down
# how many genes are significantly affected by depth
sum(res_lrt$padj < 0.1, na.rm=TRUE)
#2948
# saving them to a file
res.ordered <- res_lrt[order(res_lrt$padj),]
# adding Expression column to show the direction of change in expression, if present.
# here, the cutoff values are 0.1 for padj and 1.5 for log2FC.
res.ordered <- data.frame(res.ordered) %>%
  mutate(Expression = case_when(log2FoldChange >= log(1.5) & padj <= 0.1 ~ "Upregulated",
                                log2FoldChange <= -log(1.5) & padj <= 0.1 ~ "Downregulated",
                                TRUE ~ "Unchanged"))
head(res.ordered)
write.csv(res.ordered, file="shallow.LRT.DE.genes.csv")

# drawing a heatmap of the most significant genes
topgenes <- head(rownames(res_lrt[order(res_lrt$padj), ]), 50)
pdf("Heatmap.top.genes.shallow.pdf",width=7)
mat <- assay(vst)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds_lrt)[,c("condition","site")])
pheatmap(mat, annotation_col=df)
dev.off()

# calculating DE results for contrasts and writing them to files
# comparing S and SS
res <- results(dds, contrast=c("condition","S","SS"))
summary(res)
# 1.3 up, 8.4 down
# how many genes are DE between S and SS
sum(res$padj < 0.1, na.rm=TRUE)
# 2722
res.ordered <- res[order(res$padj),]
# adding Expression column to show the direction of change in expression, if present.
# here, the cutoff values are 0.1 for padj and 1.5 for log2FC.
res.ordered <- data.frame(res.ordered) %>%
  mutate(Expression = case_when(log2FoldChange >= log(1.5) & padj <= 0.1 ~ "Upregulated",
                           log2FoldChange <= -log(1.5) & padj <= 0.1 ~ "Downregulated",
                           TRUE ~ "Unchanged"))
head(res.ordered)
write.csv(res.ordered, file="S.SS.DE.genes.csv")

# volcano plot
#reset par
pdf("Volcano.S.SS.pdf",width=7)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot: S vs SS", xlim=c(-3,3),ylim=c(0,15)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.1 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

# comparing SS and SD
res <- results(dds, contrast=c("condition","SS","SD"))
summary(res)
# 2.1 up, 3.8 down
# how many genes are DE between SS and SD
sum(res$padj < 0.1, na.rm=TRUE)
# 1646
res.ordered <- res[order(res$padj),]
# adding Expression column to show the direction of change in expression, if present.
# here, the cutoff values are 0.1 for padj and 1.5 for log2FC.
res.ordered <- data.frame(res.ordered) %>%
  mutate(Expression = case_when(log2FoldChange >= log(1.5) & padj <= 0.1 ~ "Upregulated",
                                log2FoldChange <= -log(1.5) & padj <= 0.1 ~ "Downregulated",
                                TRUE ~ "Unchanged"))
head(res.ordered)
write.csv(res.ordered, file="SS.SD.DE.genes.csv")

# comparing S and SD
res <- results(dds, contrast=c("condition","S","SD"))
summary(res)
# 4.5 up, 15 down
# how many genes are DE between SS and SD
sum(res$padj < 0.1, na.rm=TRUE)
# 5479
res.ordered <- res[order(res$padj),]
# adding Expression column to show the direction of change in expression, if present.
# here, the cutoff values are 0.1 for padj and 1.5 for log2FC.
res.ordered <- data.frame(res.ordered) %>%
  mutate(Expression = case_when(log2FoldChange >= log(1.5) & padj <= 0.1 ~ "Upregulated",
                                log2FoldChange <= -log(1.5) & padj <= 0.1 ~ "Downregulated",
                                TRUE ~ "Unchanged"))
head(res.ordered)
write.csv(res.ordered, file="S.SD.DE.genes.csv")

# volcano plot
#reset par
pdf("Volcano.S.SD.pdf",width=7)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot: S vs SD", xlim=c(-3,3),ylim=c(0,15)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.1 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()


# Venn diagram

# writing DE expressed gene names from each comparison to a list
resSSS <- results(dds, contrast=c("condition","S","SS"))
resSSSD <- results(dds, contrast=c("condition","SS","SD"))
resSSD <- results(dds, contrast=c("condition","S","SD"))

pval_threshold <- 0.1
S.SS <- row.names(resSSS[which(resSSS$padj <= pval_threshold), ])
SS.SD <- row.names(resSSSD[which(resSSSD$padj <= pval_threshold), ])
S.SD <- row.names(resSSD[which(resSSD$padj <= pval_threshold), ])

x <- list(
  A = S.SS, 
  B = SS.SD, 
  C = S.SD)

# helper function to display Venn diagram
display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

# building Venn diagram
pdf("Venn.shallow.pdf",width=7)
display_venn(
  x,
  category.names = c("S vs SS" , "SS vs SD" , "S vs SD"),
  fill = c("#999999", "#E69F02", "#56B4E9")
)
dev.off()
#### deep ####

# reading count matrix from a file
countData  <- read.csv2('CountMatrix.csv', header=TRUE, row.names=1, sep=',')
colnames(countData) <- names
# reading metadata file
MetaData <- read.csv2('Metadata.origin.D.csv', header=TRUE, sep=',')
# removing the samples without any info in MetaData from the count table
countData <- countData[ , MetaData$id] 
# changing columns to factors
MetaData$condition <- factor(MetaData$condition, levels = c("D","DD","DS"))
MetaData$origin <- as.factor(MetaData$origin)
MetaData$batch <- as.factor(MetaData$batch)
# creating a DESeq2 object 
# since site is already nested within the origin factor, we don't need to
# additionally account for it - origin already absorbs the site effect 
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = MetaData,
                              design = ~ origin  + batch + condition)

# pre filtering - removing rows with low counts
smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# running a model
dds <- DESeq(dds)

# estimating size factors to determine if we can use vst
# to transform our data. Size factors should be less than 4 to use vst
SF <- estimateSizeFactors(dds) 
print(sizeFactors(SF))
# variance-stabilizing transformation
vst <- varianceStabilizingTransformation(dds)

# since vst does not remove variation that can be associated with covariates, 
# we manually remove the effect of covariates to be able to visualize it on PCA
mat <- assay(vst)
mm <- model.matrix(~condition, colData(vst))
mat <- limma::removeBatchEffect(mat, batch=vst$origin, batch2=vst$batch, design=mm)
assay(vst) <- mat
# PCA plot
pcaData <- plotPCA(vst, intgroup=c("condition", "site"), ntop = 1000, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("PCA.deep.pdf",width=7)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=site)) +
  geom_point(size=3) +
  ggtitle("PCA plot for translocation from deep reef") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()
dev.off()

# PCA with ellipses
pdf("PCA.deep.ellips.pdf", width=7)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  geom_mark_ellipse(aes(fill=condition), alpha=0.3) +
  ggtitle("PCA plot for translocation from deep reef") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()
dev.off()

# sample distances
sampleDists <- dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$condition, vst$site, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("Dist.deep.pdf",width=7)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

# LRT test to test the effect of depth in general
# removing "depth" to test its effect
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ origin + batch )  
res_lrt <- results(dds_lrt)
summary(res_lrt)
# 0.34 up, 2.6 down
# how many genes are significantly affected by depth
sum(res_lrt$padj < 0.1, na.rm=TRUE)
# 773
# drawing a heatmap of the most significant genes
topgenes <- head(rownames(res_lrt[order(res_lrt$padj), ]), 50)
pdf("Heatmap.top.genes.deep.pdf",width=7)
mat <- assay(vst)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds_lrt)[,c("condition","site")])
pheatmap(mat, annotation_col=df)
dev.off()
# saving them to a file
res.ordered <- res_lrt[order(res_lrt$padj),]
# adding Expression column to show the direction of change in expression, if present.
# here, the cutoff values are 0.1 for padj and 1.5 for log2FC.
res.ordered <- data.frame(res.ordered) %>%
  mutate(Expression = case_when(log2FoldChange >= log(1.5) & padj <= 0.1 ~ "Upregulated",
                                log2FoldChange <= -log(1.5) & padj <= 0.1 ~ "Downregulated",
                                TRUE ~ "Unchanged"))
head(res.ordered)
write.csv(res.ordered, file="deep.LRT.DE.genes.csv")

# calculating DE results for contrasts and writing them to files
# comparing D and DD
res <- results(dds, contrast=c("condition","D","DD"))
summary(res)
# 0.084 up, 0.061 down
# how many genes are DE between D and DD
sum(res$padj < 0.1, na.rm=TRUE)
# 38
res.ordered <- res[order(res$padj),]
# adding Expression column to show the direction of change in expression, if present.
# here, the cutoff values are 0.1 for padj and 1.5 for log2FC.
res.ordered <- data.frame(res.ordered) %>%
  mutate(Expression = case_when(log2FoldChange >= log(1.5) & padj <= 0.1 ~ "Upregulated",
                                log2FoldChange <= -log(1.5) & padj <= 0.1 ~ "Downregulated",
                                TRUE ~ "Unchanged"))
head(res.ordered)
write.csv(res.ordered, file="D.DD.DE.genes.csv")

# volcano plot
#reset par
pdf("Volcano.D.DD.pdf",width=7)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot: D vs DD", xlim=c(-3,3),ylim=c(0,15)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.1 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

# comparing DD and DS
res <- results(dds, contrast=c("condition","DD","DS"))
summary(res)
# 7% up, 0.26 down
# how many genes are DE between DD and DS
sum(res$padj < 0.1, na.rm=TRUE)
# 1917
res.ordered <- res[order(res$padj),]
# adding Expression column to show the direction of change in expression, if present.
# here, the cutoff values are 0.1 for padj and 1.5 for log2FC.
res.ordered <- data.frame(res.ordered) %>%
  mutate(Expression = case_when(log2FoldChange >= log(1.5) & padj <= 0.1 ~ "Upregulated",
                                log2FoldChange <= -log(1.5) & padj <= 0.1 ~ "Downregulated",
                                TRUE ~ "Unchanged"))
head(res.ordered)
write.csv(res.ordered, file="DD.DS.DE.genes.csv")

# volcano plot
#reset par
pdf("Volcano.DD.DS.pdf",width=7)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot: DD vs DS", xlim=c(-3,3),ylim=c(0,15)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.1 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

# comparing D and DS
res <- results(dds, contrast=c("condition","D","DS"))
summary(res)
# 12% up, 1.7 down
# how many genes are DE between DD and DS
sum(res$padj < 0.1, na.rm=TRUE)
# 3673
res.ordered <- res[order(res$padj),]
# adding Expression column to show the direction of change in expression, if present.
# here, the cutoff values are 0.1 for padj and 1.5 for log2FC.
res.ordered <- data.frame(res.ordered) %>%
  mutate(Expression = case_when(log2FoldChange >= log(1.5) & padj <= 0.1 ~ "Upregulated",
                                log2FoldChange <= -log(1.5) & padj <= 0.1 ~ "Downregulated",
                                TRUE ~ "Unchanged"))
head(res.ordered)
write.csv(res.ordered, file="D.DS.DE.genes.csv")

# volcano plot
#reset par
pdf("Volcano.D.DS.pdf",width=7)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot: D vs DS", xlim=c(-3,3),ylim=c(0,15)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.1 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

# Venn diagram

# writing DE expressed gene names from each comparison to a list
resDDD <- results(dds, contrast=c("condition","D","DD"))
resDDDS <- results(dds, contrast=c("condition","DD","DS"))
resDDS <- results(dds, contrast=c("condition","D","DS"))

pval_threshold <- 0.1
D.DD <- row.names(resSSS[which(resDDD$padj <= pval_threshold), ])
DD.DS <- row.names(resSSSD[which(resDDDS$padj <= pval_threshold), ])
D.DS <- row.names(resSSD[which(resDDS$padj <= pval_threshold), ])

x <- list(
  A = D.DD, 
  B = DD.DS, 
  C = D.DS)

# helper function to display Venn diagram
display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

# building Venn diagram
pdf("Venn.deep.pdf",width=7)
display_venn(
  x,
  category.names = c("D vs DD" , "DD vs DS" , "D vs DS"),
  fill = c("#999999", "#E69F02", "#56B4E9")
)
dev.off()


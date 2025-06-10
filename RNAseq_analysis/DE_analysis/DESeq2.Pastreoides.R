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
library(ggVennDiagram)
library(NbClust)
library(ComplexHeatmap)
library(edgeR)
library(DEGreport)
library(fgsea)
library(purrr)
library(tibble)
library(reshape2)


#load("/home/gospozha/haifa/cayman/rna/mapping/github/rin/Deseq2.rrna.rin.R")
# setting working directory 
setwd("/home/gospozha/haifa/cayman/rna/mapping/github/")

#### Preparing necessary files ####

# list files with gene counts for each sample
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

# reading count matrix from a file
countData  <- read.csv2('CountMatrix.csv', header=TRUE, row.names=1, sep=',', check.names = F)
# reading metadata file
MetaData <- read.csv2('Metadata.origin.csv', header=TRUE, sep=',')

MetaData$condition <- factor(MetaData$condition, levels = c("S", "SS", "SD", "D", "DD", "DS"))
MetaData$site <- factor(MetaData$site, levels = c("MF","CC"))
MetaData$RIN2 <- as.factor(MetaData$RIN2)
MetaData$rRNA <- as.factor(MetaData$rRNA)
MetaData$origin <- as.factor(MetaData$origin)

# testing for independence of covariates possibly related to read quality
chisq.test(MetaData$condition, MetaData$rRNA) # p value > 0.05 - they are independent
chisq.test(MetaData$condition, MetaData$RIN2) # p value > 0.05 - they are independent


#### Initial quality check ####

# Convert counts to DGEList
dge <- DGEList(counts = countData)
# remove low counts
smallestGroupSize <- 4
keep <- rowSums(dge$counts >= 10) >= smallestGroupSize
dge <- dge[keep,]

# Calculate FPM (Fragments Per Million)
fpm_values <- cpm(dge, normalized.lib.sizes = TRUE)  # edgeR's CPM is equivalent to FPM

# Convert to long format for plotting
fpm_df <- as.data.frame(fpm_values) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "FPM") %>%
  left_join(MetaData, by = c("Sample" = "id"))  # Merge with metadata

ggplot(fpm_df, aes(x = FPM, color = condition)) +
  geom_density(alpha = 0.3) +
  scale_x_log10() +
  theme_minimal() +
  labs(title="Density Plot of FPM Values per Condition",
       x="FPM (log10 scaled)")

# statistical comparison
anova_res <- aov(FPM ~ condition, data = fpm_df)
summary(anova_res)
TukeyHSD(anova_res)
# they are the same

#### DESeq2 model ####
# setting working directory 
setwd("/home/gospozha/haifa/cayman/rna/mapping/github/rin")
# creating DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = MetaData,
                              design = ~ site + rRNA + RIN2 + condition)

# pre filtering
smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds) # 31989 genes have left

# running a model
dds <- DESeq(dds)
res <- results(dds)

# Plotting histograms of p-values
hist(res$pvalue, breaks=50, col="skyblue", main="~ site + rRNA + RIN + condition",
     xlab="p-value", xlim=c(0,1), ylim=c(0, max(table(cut(res$pvalue, breaks=50)))))

# saving a DESeq2 model to an R object
saveRDS(dds, file = "dds_site_rrna_rin_condition.rds")
#dds <- readRDS(file = "../dds_site_rrna_rin_condition.rds")

#### PCA and sample distances using rlog ####

# estimating size factors to determine if it's better to use rlog
# to transform our data. rlog is more robust if size factors differ a lot
SF <- estimateSizeFactors(dds) 
print(sizeFactors(SF))

# the same using rlog transformation
rlog <- rlog(dds)
saveRDS(rlog, file = "rlog_site_rrna_rin_condition.rds")
#rlog <- readRDS(file = "rlog_site_rrna_rin_condition.rds")
# since vst does not remove variation that can be associated with covariates, 
# we manually remove the effect of covariates to be able to visualize it on PCA
mat <- assay(rlog)
mm <- model.matrix(~condition + site  + rRNA + RIN2 , colData(rlog))
treatment.design <- mm[,1:6]
batch.design <- mm[,-(1:6)]
mat <- limma::removeBatchEffect(mat, covariates=batch.design, design=treatment.design)
assay(rlog) <- mat
# PCA plot
pcaData <- plotPCA(rlog, intgroup=c("condition"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pdf("PCA.full.pdf",width=7)
pca<-ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  ggtitle("PCA of gene counts") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_hue(labels = c("10", "40",
                             "10→10","40→40",
                             "10→40", "40→10"))+
  theme_bw()
pca
ggsave("pca.jpg", pca, width = 6.5, height = 6)
norm.counts <- assay(rlog)
write.csv(norm.counts, file="./GO/rlog.counts.csv")

#dev.off()

# sample distances
sampleDists <- dist(t(assay(rlog)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rlog$condition)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#pdf("Dist.all.pdf",width=7)
dist <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dist
ggsave("dist.jpg", dist, width = 6, height = 6)
#dev.off()

#### Pairwise comparisons ####

# calculating DE results for contrasts and writing them to files
# comparing S and D
res <- results(dds, contrast=c("condition","S","D"), alpha = 0.05)
summary(res)
# 0.38 up, 0.36 down
# how many genes are DE between S and D
sum(res$padj < 0.05, na.rm=TRUE)
# 236

# adding Expression column to show the direction of change in expression, if present.
# here, the cutoff values are 0.1 for padj and 1.5 for log2FC.
res.ordered <- data.frame(res) %>%
  filter(padj<.05 & abs(log2FoldChange)>1.5)  %>%
  arrange(padj) %>%
  mutate(Expression = case_when(log2FoldChange > log(1.5) ~ "Upregulated",
                                log2FoldChange < -log(1.5) ~ "Downregulated"))
head(res.ordered)
write.csv(res.ordered, file="./GO/S.D.DE.sva1.genes.csv")

# volcano plot
#reset par
pdf("Volcano.S.D.pdf",width=7)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot: S vs D", xlim=c(-3,3),ylim=c(0,15)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

# comparing SS and DD
res <- results(dds, contrast=c("condition","SS","DD"), alpha = 0.05)
summary(res)
# 1.3 up, 0.63 down
# how many genes are DE between SS and DD
sum(res$padj < 0.05, na.rm=TRUE)
# 611
res.ordered <- data.frame(res) %>%
  filter(padj<.05 & abs(log2FoldChange)>1)  %>%
  arrange(padj) %>%
  mutate(Expression = case_when(log2FoldChange > log(1) & padj < 0.1 ~ "Upregulated",
                                log2FoldChange < -log(1) & padj < 0.1 ~ "Downregulated"))
head(res.ordered)
write.csv(res.ordered, file="./GO/SS.DD.DE.genes.lfc1.csv")

# volcano plot
#reset par
pdf("Volcano.SS.DD.pdf",width=7)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot: SS vs DD", xlim=c(-3,3),ylim=c(0,15)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

# comparing SS and DS
res <- results(dds, contrast=c("condition","SS","DS"), alpha = 0.05)
summary(res)
# 3.1 up, 0.62 down
# how many genes are DE between S and D
sum(res$padj < 0.05, na.rm=TRUE)
# 1174
res.ordered <- data.frame(res) %>%
  filter(padj<.05 & abs(log2FoldChange)>1.5)  %>%
  arrange(padj) %>%
  mutate(Expression = case_when(log2FoldChange > log(1.5) & padj < 0.1 ~ "Upregulated",
                                log2FoldChange < -log(1.5) & padj < 0.1 ~ "Downregulated"))
head(res.ordered)
write.csv(res.ordered, file="./GO/SS.DS.DE.genes.csv")

# volcano plot
#reset par
pdf("Volcano.SS.DS.pdf",width=7)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot: SS vs DS", xlim=c(-3,3),ylim=c(0,15)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

# comparing DD and SD
res <- results(dds, contrast=c("condition","DD","SD"), alpha = 0.05)
summary(res)
# 2.2 up, 4.6 down
# how many genes are DE between S and D
sum(res$padj < 0.05, na.rm=TRUE)
# 2185
res.ordered <- data.frame(res) %>%
  filter(padj<.05 & abs(log2FoldChange)>1.5)  %>%
  arrange(padj) %>%
  mutate(Expression = case_when(log2FoldChange > log(1.5) & padj < 0.1 ~ "Upregulated",
                                log2FoldChange < -log(1.5) & padj < 0.1 ~ "Downregulated"))
head(res.ordered)
write.csv(res.ordered, file="./GO/DD.SD.DE.genes.csv")

# volcano plot
#reset par
pdf("Volcano.DD.SD.pdf",width=7)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot: DD vs SD", xlim=c(-3,3),ylim=c(0,15)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

# comparing S and SS
res <- results(dds, contrast=c("condition","S","SS"), alpha = 0.05)
summary(res)
# 0.18 up, 0.58 down
# how many genes are DE
sum(res$padj < 0.05, na.rm=TRUE)
# 246
res.ordered <- data.frame(res) %>%
  filter(padj<.05 & abs(log2FoldChange)>1.5)  %>%
  arrange(padj) %>%
  mutate(Expression = case_when(log2FoldChange > log(1.5) & padj < 0.1 ~ "Upregulated",
                                log2FoldChange < -log(1.5) & padj < 0.1 ~ "Downregulated"))
head(res.ordered)
write.csv(res.ordered, file="./GO/S.SS.DE.genes.csv")

# comparing D and DD
res <- results(dds, contrast=c("condition","D","DD"), alpha = 0.05)
summary(res)
# 0.23 up, 0.58 down
# how many genes are DE
sum(res$padj < 0.05, na.rm=TRUE)
# 259
res.ordered <- data.frame(res) %>%
  filter(padj<.05 & abs(log2FoldChange)>1.5)  %>%
  arrange(padj) %>%
  mutate(Expression = case_when(log2FoldChange > log(1.5) & padj < 0.1 ~ "Upregulated",
                                log2FoldChange < -log(1.5) & padj < 0.1 ~ "Downregulated"))
head(res.ordered)
write.csv(res.ordered, file="./GO/D.DD.DE.genes.csv")


# comparing DD and DS
res <- results(dds, contrast=c("condition","DD","DS"), alpha = 0.05)
summary(res)
# 0.43 up, 0.22 down
# how many genes are DE
sum(res$padj < 0.05, na.rm=TRUE)
# 208
res.ordered <- data.frame(res) %>%
  filter(padj<.05 & abs(log2FoldChange)>1.5)  %>%
  arrange(padj) %>%
  mutate(Expression = case_when(log2FoldChange > log(1.5) & padj < 0.1 ~ "Upregulated",
                                log2FoldChange < -log(1.5) & padj < 0.1 ~ "Downregulated"))
head(res.ordered)
write.csv(res.ordered, file="./GO/DD.DS.DE.genes.csv")

# comparing S and SD
res <- results(dds, contrast=c("condition","SS","SD"), alpha = 0.05)
summary(res)
# 2.4 up, 2.2 down
# how many genes are DE
sum(res$padj < 0.05, na.rm=TRUE)
# 1487
res.ordered <- data.frame(res) %>%
  filter(padj<.05 & abs(log2FoldChange)>1.5)  %>%
  arrange(padj) %>%
  mutate(Expression = case_when(log2FoldChange > log(1.5) & padj < 0.1 ~ "Upregulated",
                                log2FoldChange < -log(1.5) & padj < 0.1 ~ "Downregulated"))
head(res.ordered)
write.csv(res.ordered, file="./GO/SS.SD.DE.genes.csv")


#### Venn diagrams ####

res7 <- results(dds, contrast=c("condition","S","SS"), alpha = 0.05)
summary(res7)
# up 0.18%, down 0.32%
res8 <- results(dds, contrast=c("condition","D","DD"), alpha = 0.05)
summary(res8)
# up 11%, down 0.9
res1 <- results(dds, contrast=c("condition","S","D"), alpha = 0.05)
summary(res1)
# up 0.25%, down 0.52%
res2 <- results(dds, contrast=c("condition","SS","DD"), alpha = 0.05)
summary(res2)
# heat stress related genes?
# up 9.7%, down 1.5%
res3 <- results(dds, contrast=c("condition","DD","SD"), alpha = 0.05)
summary(res3)
# up 0.42%, down 1.7%
res4 <- results(dds, contrast=c("condition","SS","DS"), alpha = 0.05)
summary(res4)
# up 11%, down 1.9%
res5 <- results(dds, contrast=c("condition","DD","DS"), alpha = 0.05)
summary(res5)
# up 0.07, down 0,01%
res6 <- results(dds, contrast=c("condition","SS","SD"), alpha = 0.05)
summary(res6)
# up 2.6, down 0.79
res9 <- results(dds, contrast=c("condition","SD","DS"), alpha = 0.05)
summary(res9)
# up 3.2%, down 0.75
res10 <- results(dds, contrast=c("site","CC","MF"), alpha = 0.05)
summary(res10)
# up 0%, down 0

pval_threshold <- 0.05
S.D <- row.names(res1[which(res1$padj <= pval_threshold), ])
SS.DD <- row.names(res2[which(res2$padj <= pval_threshold), ])
DD.SD <- row.names(res3[which(res3$padj <= pval_threshold), ])
SS.DS <- row.names(res4[which(res4$padj <= pval_threshold), ])
DD.DS <- row.names(res5[which(res5$padj <= pval_threshold), ])
SS.SD <- row.names(res6[which(res6$padj <= pval_threshold), ])
S.SS <- row.names(res7[which(res7$padj <= pval_threshold), ])
D.DD <- row.names(res8[which(res8$padj <= pval_threshold), ])
SD.DS <- row.names(res9[which(res9$padj <= pval_threshold), ])


x <- list(
  A = S.D, 
  B = SS.DD, 
  G = S.SS,
  H = D.DD)

cont <- ggVennDiagram(x,  
                      category.names = c("S vs D", 'SS vs DD', 'S vs SS', 'D vs DD'),
                      label_alpha = 0, label = "count")+
  scale_fill_distiller(palette = "OrRd", direction = 1)
cont
ggsave("./GO/venn.contr.jpg", cont, width = 7, height = 6)

x <- list(
  A = SS.DD, 
  B = SS.DS, 
  G = SS.SD,
  H = DD.SD,
  I = DD.DS)

treat <- ggVennDiagram(x,  
                       category.names = c("SS vs DD", 'SS vs DS', 'SS vs SD', 
                                          'DD vs SD', 'DD vs DS'),
                       label_alpha = 0, label = "count")+
  scale_fill_distiller(palette = "OrRd", direction = 1)
treat
ggsave("./GO/venn.treat.jpg", treat, width = 7, height = 6)

x <- list(
  A = SS.DD, 
  B = SS.DS, 
  С = DD.SD)

treat <- ggVennDiagram(x,  
                       category.names = c("SS vs DD", 'SS vs DS', "DD vs SD"),
                       label_alpha = 0, label = "count")+
  scale_fill_distiller(palette = "OrRd", direction = 1)
treat
ggsave("./GO/venn.DDcompare.jpg", treat, width = 7, height = 6)


x <- list(
  A = SS.SD, 
  B = SS.DS, 
  С = DD.DS,
  D = DD.SD,
  E = SD.DS)

treat <- ggVennDiagram(x,  
                       category.names = c("SS vs SD", 'SS vs DS', "DD vs DS",
                                          "DD vs SD", "SD vs DS"),
                       label_alpha = 0, label = "count")+
  scale_fill_distiller(palette = "OrRd", direction = 1)
treat
ggsave("./GO/venn.trans.jpg", treat, width = 7, height = 6)

#### visualizing the results ####
contrast_list <- list(
  c("condition", "S", "D"),
  c("condition", "D", "DD"),
  c("condition", "S", "SS"),
  c("condition", "SS", "DD"),
  c("condition", "SS", "SD"),
  c("condition", "SS", "DS"),
  c("condition", "DD", "SD"),
  c("condition", "DD", "DS")
  #c("condition", "SD", "DS")
)
# Placeholder list to store results
deg_summary_list <- list()

for (con in contrast_list) {
  res <- results(dds, contrast = con, alpha = 0.05)
  
  # Remove NAs
  res <- na.omit(res)
  
  # Count
  up <- sum(res$padj < 0.05 & res$log2FoldChange > 0)
  down <- sum(res$padj < 0.05 & res$log2FoldChange < 0)
  total <- up + down
  
  deg_summary_list[[paste(con[2], "vs", con[3])]] <- data.frame(
    Contrast = paste(con[2], "vs", con[3]),
    Up = up,
    Down = down,
    TotalDEGs = total
  )
}

# Combine into single data frame
df_summary <- bind_rows(deg_summary_list)

df_long <- df_summary %>%
  tidyr::pivot_longer(cols = c("Up", "Down"), names_to = "Direction", values_to = "Count") %>%
  mutate(Percent = round(100 * Count / 31989, 1))
df_long$Direction <- factor(df_long$Direction, levels = c("Up", "Down"))

deg.bar <- ggplot(df_long, aes(x = Contrast, y = Count, fill = Direction)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Percent, "%")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black") +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.1, 0.1))) +
  scale_fill_manual(values = c("Down" = "aquamarine3", "Up" = "coral")) +
  theme_bw() +
  labs(y = "Number of DEGs",
    x = "Contrast",
    fill = "Direction",
    title = "Up- and downregulated genes per contrast",
    subtitle = "p.adj<0.05") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
deg.bar
ggsave("./GO/degs.bar.jpg", deg.bar, width = 6, height = 6)


#### LRT ####
# LRT test to test the effect of depth in general
# removing "depth" to test its effect
dds_lrt <- DESeq(dds, test="LRT", reduced =  ~ site + rRNA + RIN2)  
saveRDS(dds_lrt, "depth.model.dds.R")
res_lrt <- results(dds_lrt, alpha = 0.05)
summary(res_lrt)
# 5.9 up, 4.8 down
# how many genes are significantly affected by depth
sum(res_lrt$padj < 0.05, na.rm=TRUE)
# 3414
# saving them to a file
res.ordered <- data.frame(res) %>%
  filter(padj<.05)  %>%
  arrange(padj)

head(res.ordered)
write.csv(res.ordered, file="./GO/all.cond.DE.genes.csv")

# drawing a heatmap of the most significant genes
topgenes <- head(rownames(res_lrt[order(res_lrt$padj), ]), 200)
mat <- assay(rlog)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds_lrt)[,c("condition")])
rownames(df) <- colnames(mat)
pheatmap(mat, annotation_col=df)

topgenes <- rownames(res_lrt[order(res_lrt$padj), ])
mat <- assay(rlog)[topgenes,]
mat <- mat - rowMeans(mat)
# Find best number of clusterw
set.seed(1990)
# list of potential indices (note: I removed some which took longer to run)
indices <- c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "dunn", "sdindex", "sdbw")
# initialize var to collect clustering results
results <- list()
# loop over indices with try function (which continues running even with errors)
for (i in 1:length(indices)) {
  print(paste0("Trying ", indices[i], " index..."))
  results[[i]] <- try(NbClust(data=mat,min.nc=2,max.nc=15, index=indices[i], method="kmeans")) 
}

num_clust <- list()
for (i in 1:length(results)){
  num_clust[[i]] <- try(as.numeric(results[[i]]$Best.nc[1]))
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

paste0("Based on a number of criteria, we will select ", getmode(num_clust), " clusters.")
# 2
calc.kmeans <- kmeans(mat, 2)
cluster_res <- data.frame(gene_id = names(calc.kmeans$cluster), cluster = calc.kmeans$cluster)
mat2 <- as.data.frame(assay(rlog)[topgenes,])
mat2$gene_id <- rownames(mat2)
pln_DEGs_all_clust <- merge(mat2, cluster_res, by = "gene_id")

pln_DEGs_all_clust <- subset(pln_DEGs_all_clust, select = c(gene_id, cluster))

#Prepare annotations
hm_ann_row <- unique(pln_DEGs_all_clust)
rownames(hm_ann_row) <- hm_ann_row$gene_id
hm_ann_row <- subset(hm_ann_row, select=cluster)
hm_ann_row$cluster <- gsub(1,"Cluster1",hm_ann_row$cluster)
hm_ann_row$cluster <- gsub(2,"Cluster2",hm_ann_row$cluster)
hm_ann_row <- as.matrix(hm_ann_row[rownames(mat),])

hmTreatment <- colData(rlog)[c("condition")]

hm_ann_col <- ComplexHeatmap::HeatmapAnnotation(df=hmTreatment, 
                                col = list(treatment=c("D" ="blue", 
                                                       "S" ="green",
                                                       "DD"  ="red",
                                                       "SS" = "grey",
                                                       "DS" = "pink",
                                                       "SD" = "aquamarine"))) #make dataframe for column naming

pln_DEGheatmap <-  Heatmap(mat, column_title = "Treatment", 
                           name = "expression",
                           show_row_names = FALSE, top_annotation = hm_ann_col, show_column_names = FALSE, row_dend_side = "left" ,
                           column_split = 6, column_dend_height = unit(0.5, "in"),
                           km = 2, row_km_repeats = 100, row_title = c("Cluster1", "Cluster2"),
                           row_gap = unit(2.5, "mm"), border = TRUE,
                           column_names_gp =  gpar(fontsize = 10)); pln_DEGheatmap

pln_DEGheatmap
ggsave("./GO/depth.heatmap.jpg", pln_DEGheatmap, width = 12, height = 6)

# degPatterns

# Obtain rlog values for those significant genes
cluster_rlog <- assay(rlog)[rownames(res.ordered), ]

rownames(MetaData) <- MetaData$id
# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
clusters <- DEGreport::degPatterns(cluster_rlog, metadata = MetaData, 
                                   time = "condition")
# Extract the Group 1 genes
cluster_groups <- clusters$df
group1 <- clusters$df %>%
  filter(cluster == 1)

save.image("Deseq2.rrna.rin.R")


#### biomineralization-related genes ####
setwd("/home/gospozha/haifa/cayman/rna/biomin/")

# heatmap of Z-scores
biomineralization_gene_list <- read.csv("past.biomin.genes.txt", sep = " ", header = F) 
bio_genes <- c(biomineralization_gene_list$V1) 
bio_genes <- bio_genes[bio_genes %in% rownames(rlog)]
mat <- assay(rlog)[bio_genes, ]
mat_z <- t(scale(t(mat)))  # Z-score by gene
anno_col <- data.frame(condition = colData(rlog)$condition)
rownames(anno_col) <- colnames(rlog)  

# heatmap of selected genes from rlog
pheatmap(mat_z,
         annotation_col = anno_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE)  # optional for clean plots


# are these genes significantly participate in depth change?
res.ordered <- data.frame(res_lrt) %>%
  filter(padj<.05)  %>%
  arrange(padj)
lrt.biomin <- res.ordered[bio_genes, ]
lrt.biomin <- na.omit(lrt.biomin)
write.csv(lrt.biomin, file="./GO/lrt.biomin.genes.csv")

# visualizing the boxplots for these genes
sig_bio_genes <- rownames(lrt.biomin) 

# Make sure they exist in rld
sig_bio_genes <- sig_bio_genes[sig_bio_genes %in% rownames(rlog)]

# Extract rlog expression matrix for those genes
expr_mat <- assay(rlog)[sig_bio_genes, ]

# Transpose and convert to data.frame
df <- as.data.frame(t(expr_mat))
df$sample <- rownames(df)
df$condition <- colData(rlog)$condition[match(df$sample, rownames(colData(rlog)))]

# Pivot longer for ggplot
df_long <- df %>%
  pivot_longer(cols = all_of(sig_bio_genes),
               names_to = "gene",
               values_to = "expression")

# Plot
ggplot(df_long, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(title = "rlog Expression of biomineralization genes significant for depth",
       y = "rlog Expression",
       x = "Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


## GSEA analysis

# 1. Create a named list of your DESeq2 results
res_list <- list(
  S_D = res1,
  SS_DD = res2,
  DD_SD = res3,
  SS_DS = res4,
  DD_DS = res5,
  SS_SD = res6,
  S_SS = res7,
  D_DD = res8
)

# 2. Your biomineralization gene set (named list)
bio_gene_set <- list("Biomineralization" = bio_genes)

# 3. Function to run fgsea on a DESeq2 result
run_fgsea <- function(res, gene_set, name) {
  res_df <- as.data.frame(res)
  
  # Remove rows with NA values in LFC or padj
  res_df <- res_df %>%
    filter(!is.na(log2FoldChange), !is.na(padj))
  
  # Compute ranking score - rank according to pval and logfc
  res_df$rank_score <- res_df$log2FoldChange * -log10(res_df$padj + 1e-300) 
  
  # Create named vector
  ranked <- res_df$rank_score
  names(ranked) <- rownames(res_df)
  
  # Remove infinite or NaN values (can happen if padj = 0)
  ranked <- ranked[is.finite(ranked)]
  
  # Sort in decreasing order
  ranked <- sort(ranked, decreasing = TRUE)
  
  # Run fgsea
  fgsea_res <- fgsea(pathways = gene_set, stats = ranked)
  fgsea_res$contrast <- name
  fgsea_res
}


# 4. Run fgsea on all DE results and bind results into one table
all_fgsea <- imap_dfr(res_list, ~ run_fgsea(.x, bio_gene_set, .y))

# 5. Filter significant results (FDR < 0.05)
signif_fgsea <- all_fgsea %>%
  filter(padj < 0.05)

# 6. Optional: format for presentation
signif_fgsea %>%
  select(contrast, pathway, NES, padj, leadingEdge, size) %>%
  arrange(contrast)

## pathway is not significantly enriched in any of contrasts

## Are these genes significantly DE in the contrasts?

# Function to extract direction of significant biomineralization genes from a single res
get_directional_sig_bio_genes <- function(res, bio_genes) {
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Filter biomineralization genes with padj < 0.05
  res_df <- res_df %>%
    filter(gene %in% bio_genes, !is.na(padj), padj < 0.05)
  
  # Determine direction of change
  res_df <- res_df %>%
    mutate(direction = case_when(
      log2FoldChange > 0 ~ "Up",
      log2FoldChange < 0 ~ "Down",
      TRUE ~ "0"
    )) %>%
    select(gene, direction)
  
  return(res_df)
}

# Apply to all results, get a named list of data.frames
direction_lists <- lapply(res_list, get_directional_sig_bio_genes, bio_genes = bio_genes)
names(direction_lists) <- names(res_list)

# Get all genes that were significant in at least one contrast
all_sig_genes <- unique(unlist(lapply(direction_lists, \(df) df$gene)))

# Build a gene × contrast matrix filled with "0"
summary_df <- matrix("0", nrow = length(all_sig_genes), ncol = length(res_list),
                     dimnames = list(all_sig_genes, names(res_list)))

# Fill in "Up" or "Down" for significant cases
for (contrast_name in names(direction_lists)) {
  df <- direction_lists[[contrast_name]]
  summary_df[df$gene, contrast_name] <- df$direction
}

# Convert to data frame with gene column first
summary_df <- as.data.frame(summary_df)
summary_df$gene <- rownames(summary_df)
summary_df <- summary_df %>% select(gene, everything())

# View it
print(summary_df)


write.csv(summary_df, "/home/gospozha/haifa/cayman/rna/biomin/biomineralization_gene_presence_summary.csv", row.names = FALSE)

## draw them
# Suppose your significant genes vector is called:
sig_bio_genes <- summary_df$gene  # or your specific 14 genes

# Make sure they exist in rld
sig_bio_genes <- sig_bio_genes[sig_bio_genes %in% rownames(rlog)]

# Extract rlog expression matrix for those genes
expr_mat <- assay(rlog)[sig_bio_genes, ]

# Transpose and convert to data.frame
df <- as.data.frame(t(expr_mat))
df$sample <- rownames(df)
df$condition <- colData(rlog)$condition[match(df$sample, rownames(colData(rlog)))]

# Pivot longer for ggplot
df_long <- df %>%
  pivot_longer(cols = all_of(sig_bio_genes),
               names_to = "gene",
               values_to = "expression")

# Plot
ggplot(df_long, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(title = "rlog Expression of Significant Biomineralization Genes across contrasts",
       y = "rlog Expression",
       x = "Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

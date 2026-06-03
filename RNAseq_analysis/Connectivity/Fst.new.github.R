library(gdsfmt)
library(SNPRelate)
library(hierfstat)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)
library(gridGraphics)
library(cowplot)
library(StAMPP)
library(adegenet)
library(ggrepel)
library(GWASTools)
library(poppr)
library(patchwork)
library(ape)
library(ggtree)
library(pcadapt)
library(qvalue)

setwd("~/haifa/cayman/rna/connect/admixture/new/revision")


#### no artificial clones ####
# all except for the artificial clones and samples below missing 0.9 threshold
## Load the pruned SNPs

bed.fn <- "P_astreoides_vcf_pruned_no_art_0.4.bed"
fam.fn <- "P_astreoides_vcf_pruned_no_art_0.4.fam"
bim.fn <- "P_astreoides_vcf_pruned_no_art_0.4.bim"

## Create gds file
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "plink_pruned_no_art.gds", cvt.chr="char")
snpgdsSummary("plink_pruned_no_art.gds")
## Open the gds file
genofile <- snpgdsOpen("plink_pruned_no_art.gds")
#The total number of samples: 20 
#The total number of SNPs: 6337 
###

## Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Load metadata 
meta <- read.csv("Metadata.csv")
# Ensure sample names match
meta <- meta[match(sample.id, meta$id), ]

geno <- snpgdsGetGeno(genofile)  # matrix of SNPs x samples

# hierfstat expects genotypes in 1/2 digit allele coding (e.g., 11, 12, 22), not 0/1/2

#Convert 0/1/2 to hierfstat genotype format
geno.df <- as.data.frame(t(geno))

# Convert 0/1/2 to hierfstat genotype format
geno.df.hierf <- apply(geno.df, 2, function(x) {
  y <- rep(NA, length(x))
  y[x == 0] <- 11
  y[x == 1] <- 12
  y[x == 2] <- 22
  return(y)
})
geno.df.hierfstat <- as.data.frame(t(geno.df.hierf))


geno.df.hierfstat$pop <- meta$without.relatives
# hierfstat wants pop column to be first
geno.df.hierfstat <- geno.df.hierfstat[, c(ncol(geno.df.hierfstat), 1:(ncol(geno.df.hierfstat) - 1))]

genotype <- geno.df.hierfstat
genotype.AA <- apply(t(genotype), 2, function(x) {
  y <- rep(NA, length(x))
  y[x == 11] <- ("AA")
  y[x == 12] <- ("AB")
  y[x == 22] <- ("BB")
  return(y)
})
genotype <- as.data.frame(t(genotype.AA))

genotype$Sample <- as.factor(sample.id)
genotype$Ploidy <-2
genotype$Format <- "BiA"
genotype$Pop <- meta$without.relatives
genotype <- genotype[!is.na(genotype$Pop), ]

genotype <- genotype[,c("Sample", "Pop", "Ploidy", "Format", setdiff(names(genotype), c("Sample", "Pop", "Ploidy", "Format") ))]
cols <- 4:48477
genotype[] <-lapply(genotype, as.factor)
genotype$Ploidy <- 2
genotype.st <- stamppConvert(genotype, "r")


#### PCA ####
genotype.st  -> genotype.st2
geno.nas <- genotype.st2[, !(names(genotype.st2) %in% c("Sample", "Pop", "pop.num", "ploidy", "format"))]
pop.vector <- as.factor(genotype.st2[[2]])
convert_to_genind_format <- function(x) {
  ifelse(is.na(x), NA,
         ifelse(x == 1, "11",
                ifelse(x == 0.5, "12",
                       ifelse(x == 0, "22", NA))))
}

geno.char <- as.data.frame(lapply(geno.nas, convert_to_genind_format), stringsAsFactors = FALSE)
geno.obj <- df2genind(geno.char, pop = pop.vector, ploidy = 2, NA.char = NA, sep = "")
dist.mtx <- dist(geno.obj)
pca <- dudi.pca(tab(geno.obj, NA.method = "mean"), scannf = F, nf = 2)
s.class(pca$li, fac = pop.vector, col = rainbow(length(unique(pop.vector))))


# Extract individual scores 
scores <- as.data.frame(pca$li)

# Ensure scores$SampleID contains names like "./35-MF-SD" instead of "1", "2", "3"
scores$SampleID <- genotype.st2$Sample 

scores$Pop <- pop.vector

# Variance explained
eig <- 100 * pca$eig / sum(pca$eig)

# Calculate centroids per population
centroids <- scores %>%
  group_by(Pop) %>%
  summarise(
    Axis1 = mean(Axis1),
    Axis2 = mean(Axis2)
  )

# Add centroid coords to scores
scores <- scores %>%
  left_join(centroids, by = "Pop", suffix = c("", ".centroid"))

# Shorten names to numbers only
scores$label <- gsub("^\\./|\\-.*$", "", scores$SampleID)

# inch to pt conversion 
pt_to_mm <- function(pt) pt / 2.845

# Plot
pca_plot_no_art <- ggplot(scores, aes(x = Axis1, y = Axis2, color = Pop)) +
  # segments from sample to centroid
  geom_segment(aes(xend = Axis1.centroid, yend = Axis2.centroid),
               alpha = 0.4, linewidth = 0.3) +
  
  # centroid labels
  geom_text(data = centroids, aes(x = Axis1, y = Axis2, label = Pop, color = Pop),
            size = pt_to_mm(8), show.legend = FALSE, alpha = 0.6, fontface = "bold") +
  
  # points for samples
  geom_point(size = 1.5) +
  
  # sample numbers
  geom_text_repel(aes(label = label), 
                  size = pt_to_mm(6), 
                  segment.size = 0.28,
                  segment.alpha = 0.6,
                  show.legend = FALSE,
                  max.overlaps = Inf) + 
  
  labs(
    x = paste0("PC1 (", round(eig[1], 1), "%)"),
    y = paste0("PC2 (", round(eig[2], 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "none")


pca_plot_no_art
saveRDS(pca_plot_no_art, "fig_pca_no_art.rds")
pca_plot_no_art <- readRDS("fig_pca_no_art.rds")

#### IBS dendrogram and tree ####

# Identity by State (IBS) Analysis
ibs <- snpgdsIBS(genofile, num.thread=2)
set.seed(100)
ibs.hc <- snpgdsHCluster(ibs)
rv <- snpgdsCutTree(ibs.hc)

#  IBS similarity matrix into a distance matrix
dist_matrix <- as.dist(1 - ibs$ibs)

# hierarchical clustering (UPGMA)
standard_hclust <- hclust(dist_matrix, method = "average")

# sample IDs to the cluster labels
standard_hclust$labels <- ibs$sample.id

#  phylo format
ibs_phylo <- as.phylo(standard_hclust)

# Shorten the tip labels in the phylo object to numbers only
# Removes "./" and everything after the first "-"
ibs_phylo$tip.label <- gsub("^\\./|\\-.*$", "", ibs_phylo$tip.label)

# Prepare metadata for ggtree
meta_subset <- meta[, c("id", "all")] 
colnames(meta_subset)[1] <- "label"

# Ensure metadata labels also match the format
meta_subset$label <- gsub("^\\./|\\-.*$", "", meta_subset$label)

# 5. Plot with ggtree
ibs_no_art <- ggtree(ibs_phylo, layout="rectangular", linewidth = 0.3) %<+% meta_subset +
  geom_tiplab(size=pt_to_mm(6), offset=0.01) + 
  geom_tippoint(aes(color=all), size=1.5) +
  theme_tree2() +
  labs(color = "Group")+
  xlim(0, 0.21)

ibs_no_art
saveRDS(ibs_no_art, "fig_ibs_no_art.rds")
ibs_no_art <- readRDS("fig_ibs_no_art.rds")

#### no relatives ####
# all clones and first-degree relatives detected by KING distance removed
## Load the pruned SNPs

bed.fn <- "P_astreoides_vcf_pruned_no_rel_0.4.bed"
fam.fn <- "P_astreoides_vcf_pruned_no_rel_0.4.fam"
bim.fn <- "P_astreoides_vcf_pruned_no_rel_0.4.bim"

## Create gds file
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "plink_pruned_no_rel.gds", cvt.chr="char")
snpgdsSummary("plink_pruned_no_rel.gds")
## Open the gds file
genofile <- snpgdsOpen("plink_pruned_no_rel.gds")
#The total number of samples: 13 
#The total number of SNPs: 6337 
###

## Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Load metadata 
meta <- read.csv("Metadata.csv")
# Ensure sample names match
meta <- meta[match(sample.id, meta$id), ]

geno <- snpgdsGetGeno(genofile)  # matrix of SNPs x samples

#Convert 0/1/2 to hierfstat genotype format
geno.df <- as.data.frame(t(geno))

# Convert 0/1/2 to hierfstat genotype format
geno.df.hierf <- apply(geno.df, 2, function(x) {
  y <- rep(NA, length(x))
  y[x == 0] <- 11
  y[x == 1] <- 12
  y[x == 2] <- 22
  return(y)
})
geno.df.hierfstat <- as.data.frame(t(geno.df.hierf))


geno.df.hierfstat$pop <- meta$without.relatives
# hierfstat wants pop column to be first
geno.df.hierfstat <- geno.df.hierfstat[, c(ncol(geno.df.hierfstat), 1:(ncol(geno.df.hierfstat) - 1))]

genotype <- geno.df.hierfstat
genotype.AA <- apply(t(genotype), 2, function(x) {
  y <- rep(NA, length(x))
  y[x == 11] <- ("AA")
  y[x == 12] <- ("AB")
  y[x == 22] <- ("BB")
  return(y)
})
genotype <- as.data.frame(t(genotype.AA))

genotype$Sample <- as.factor(sample.id)
genotype$Ploidy <-2
genotype$Format <- "BiA"
genotype$Pop <- meta$without.relatives
genotype <- genotype[!is.na(genotype$Pop), ]

genotype <- genotype[,c("Sample", "Pop", "Ploidy", "Format", setdiff(names(genotype), c("Sample", "Pop", "Ploidy", "Format") ))]
cols <- 4:48477
genotype[] <-lapply(genotype, as.factor)
genotype$Ploidy <- 2
genotype.st <- stamppConvert(genotype, "r")

#### Fst ####

genotype.fst <- stamppFst(genotype.st,  nboots = 1000, percent = 95, nclusters=3)
genotype.fst[["Pvalues"]]
genotype.fst[["Fsts"]]

# confidence intervals
# Extract the bootstrap values (columns 3 onwards)
boot_data <- genotype.fst[["Bootstraps"]][, -c(1, 2)]

# Calculate the 2.5% and 97.5% quantiles for each row
ci_lower <- apply(boot_data, 1, quantile, probs = 0.025, na.rm = TRUE)
ci_upper <- apply(boot_data, 1, quantile, probs = 0.975, na.rm = TRUE)

# Combine into a clean summary table
fst_summary <- data.frame(
  Pop1 = genotype.fst[["Bootstraps"]]$Population1,
  Pop2 = genotype.fst[["Bootstraps"]]$Population2,
  Fst = genotype.fst[["Fsts"]][lower.tri(genotype.fst[["Fsts"]])], # Extracts lower triangle values
  L95 = ci_lower,
  U95 = ci_upper,
  Pval = genotype.fst[["Pvalues"]][lower.tri(genotype.fst[["Pvalues"]])]
)

print(fst_summary)

#### plot Fst no relatives ####

# Prepare data
fst_mat  <- genotype.fst[["Fsts"]]
# Ensure p-values are available for both triangles by filling the matrix
pval_mat <- genotype.fst[["Pvalues"]]
boot_data <- genotype.fst[["Bootstraps"]][, -c(1, 2)]

# Calculate CIs
ci_low <- apply(boot_data, 1, quantile, probs = 0.025, na.rm = TRUE)
ci_upp <- apply(boot_data, 1, quantile, probs = 0.975, na.rm = TRUE)

# Convert to Long Format
fst_df <- expand.grid(Pop1 = rownames(fst_mat), 
                      Pop2 = colnames(fst_mat), 
                      stringsAsFactors = FALSE)

# Add values and labels
fst_df$label   <- ""
fst_df$color_val <- NA_real_
boot_idx <- 1

for(i in 1:nrow(fst_df)) {
  p1 <- fst_df$Pop1[i]
  p2 <- fst_df$Pop2[i]
  
  r <- which(rownames(fst_mat) == p1)
  c <- which(colnames(fst_mat) == p2)
  
  if (r > c) {
    # LOWER TRIANGLE: Fst [CI]
    fst_val <- fst_mat[r, c]
    fst_df$label[i] <- sprintf("%.2f\n[%.2f, %.2f]", fst_val, ci_low[boot_idx], ci_upp[boot_idx])
    fst_df$color_val[i] <- fst_val
    boot_idx <- boot_idx + 1
  } else if (r < c) {
    # UPPER TRIANGLE: P-values 
    pv <- pval_mat[c, r] 
    fst_df$label[i] <- ifelse(is.na(pv), "", 
                              ifelse(pv < 0.001, "p < 0.001", sprintf("p = %.3f", pv)))
  }
}


#  Final plot 
fst_no_rel <- ggplot(fst_df, aes(x = Pop2, y = Pop1)) +
  geom_tile(fill = "white", color = "grey90", lwd = 0.3) +
  geom_tile(aes(fill = color_val), color = "grey90", lwd = 0.3) +
  
  geom_text(
    aes(label = label),
    size = pt_to_mm(6),
    color = "black",
    lineheight = 0.75
  ) +
  
  scale_fill_gradientn(
    colors = c("#e0f3f8", "#abd9e9", "#74add1", "#4575b4"),
    name = expression(F[ST]),
    na.value = "transparent",
    guide = guide_colorbar(
      barwidth = unit(3, "mm"),
      barheight = unit(22, "mm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  
  scale_y_discrete(limits = rev(rownames(fst_mat))) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    
    axis.text = element_text(
      face = "bold",
      color = "black",
      size = 6
    ),
    
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.width = unit(3, "mm"),
    legend.key.height = unit(4, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, -4),
    
    plot.margin = margin(3, 1, 3, 3, unit = "pt")
  )
fst_no_rel
saveRDS(fst_no_rel, "fig_fst_no_rel.rds")
fst_no_rel <- readRDS("fig_fst_no_rel.rds")

#### PCA_no_rel ####
genotype.st  -> genotype.st2
geno.nas <- genotype.st2[, !(names(genotype.st2) %in% c("Sample", "Pop", "pop.num", "ploidy", "format"))]
pop.vector <- as.factor(genotype.st2[[2]])
convert_to_genind_format <- function(x) {
  ifelse(is.na(x), NA,
         ifelse(x == 1, "11",
                ifelse(x == 0.5, "12",
                       ifelse(x == 0, "22", NA))))
}

geno.char <- as.data.frame(lapply(geno.nas, convert_to_genind_format), stringsAsFactors = FALSE)
geno.obj <- df2genind(geno.char, pop = pop.vector, ploidy = 2, NA.char = NA, sep = "")
dist.mtx <- dist(geno.obj)
pca <- dudi.pca(tab(geno.obj, NA.method = "mean"), scannf = F, nf = 2)
s.class(pca$li, fac = pop.vector, col = rainbow(length(unique(pop.vector))))

scores <- as.data.frame(pca$li)

scores$SampleID <- genotype.st2$Sample 

scores$Pop <- pop.vector

eig <- 100 * pca$eig / sum(pca$eig)

# centroids per population
centroids <- scores %>%
  group_by(Pop) %>%
  summarise(
    Axis1 = mean(Axis1),
    Axis2 = mean(Axis2)
  )

# centroid coords to scores
scores <- scores %>%
  left_join(centroids, by = "Pop", suffix = c("", ".centroid"))

# truncate names
scores$label <- gsub("^\\./|\\-.*$", "", scores$SampleID)

# 2. Plot
pca_plot_no_rel <- ggplot(scores, aes(x = Axis1, y = Axis2, color = Pop)) +
  # segments from sample to centroid
  geom_segment(aes(xend = Axis1.centroid, yend = Axis2.centroid),
               alpha = 0.4, linewidth = 0.3) +
  
  # centroid labels
  geom_text(data = centroids, aes(x = Axis1, y = Axis2, label = Pop, color = Pop),
            size = pt_to_mm(8), show.legend = FALSE, alpha = 0.6, fontface = "bold") +
  
  # points for samples
  geom_point(size = 1.5) +
  
  # Sample numbers
  geom_text_repel(aes(label = label), 
                  size = pt_to_mm(6), 
                  segment.size = 0.28,
                  segment.alpha = 0.6,
                  show.legend = FALSE,
                  max.overlaps = Inf) + 
  
  labs(
    x = paste0("PC1 (", round(eig[1], 1), "%)"),
    y = paste0("PC2 (", round(eig[2], 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "none")

pca_plot_no_rel


pca_plot_no_rel
saveRDS(pca_plot_no_rel, "fig_pca_no_rel.rds")
pca_plot_no_rel <- readRDS("fig_pca_no_rel.rds")


#### IBS dendrogram and tree no relatives ####
# Identity by State (IBS) Analysis
ibs <- snpgdsIBS(genofile, num.thread=2)
set.seed(100)
ibs.hc <- snpgdsHCluster(ibs)
rv <- snpgdsCutTree(ibs.hc)
dist_matrix <- as.dist(1 - ibs$ibs)
standard_hclust <- hclust(dist_matrix, method = "average")
standard_hclust$labels <- ibs$sample.id
ibs_phylo <- as.phylo(standard_hclust)
ibs_phylo$tip.label <- gsub("^\\./|\\-.*$", "", ibs_phylo$tip.label)
meta_subset <- meta[, c("id", "all")] 
colnames(meta_subset)[1] <- "label"
meta_subset$label <- gsub("^\\./|\\-.*$", "", meta_subset$label)

ibs_no_rel <-ggtree(ibs_phylo, layout="rectangular", linewidth = 0.3) %<+% meta_subset +
  geom_tiplab(size=pt_to_mm(6), offset=0.01) + 
  geom_tippoint(aes(color=all), size=1.5) +
  theme_tree2() +
  labs(color = "Group")+
  xlim(0, 0.21)

ibs_no_rel
saveRDS(ibs_no_rel, "fig_ibs_no_rel.rds")
ibs_no_rel <- readRDS("fig_ibs_no_rel.rds")

#### pcadapt ####

# filtered plink files 
path_to_bed <- "P_astreoides_vcf_pruned_no_rel_0.4.bed"

# load the data
# pcadapt will look for the associated .bim and .fam files in the same folder
data_pcadapt <- read.pcadapt(path_to_bed, type = "bed")

# analysis
x <- pcadapt(data_pcadapt, K = 11)
plot(x, option = "screeplot")
# K = 3 
x <- pcadapt(data_pcadapt, K = 3)

plot(x, option = "scores",  i = 2, j = 3, pop = meta$all)
# plots
plot(x, option = "scores")
plot(x, option = "qqplot")
plot(x, option = "stat.distribution")
plot(x, option = "manhattan")

# checking different K
Ks <- 1:4

res_list <- lapply(Ks, function(k) {
  x <- pcadapt(data_pcadapt, K = k, ploidy = 2)
  q <- qvalue::qvalue(x$pvalues)$qvalues
  data.frame(
    K = k,
    n_outliers_q01 = sum(q < 0.01, na.rm = TRUE),
    n_outliers_q05 = sum(q < 0.05, na.rm = TRUE)
  )
})

do.call(rbind, res_list)

# significance threshold
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05 
outliers <- which(qval < alpha)

# Extract names of outlier SNPs from  .bim file
bim <- read.table("P_astreoides_vcf_pruned_no_rel_0.4.bim")
outlier_snp_names <- bim$V2[outliers]
# 529
# Save for PLINK
write.table(outlier_snp_names, "outliers_to_remove_0.05.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

#### no relatives, pcadapt snps removed ####

## Load the pruned SNPs

bed.fn <- "P_astreoides_vcf_pruned_no_rel_0.4_pcadapt.bed"
fam.fn <- "P_astreoides_vcf_pruned_no_rel_0.4_pcadapt.fam"
bim.fn <- "P_astreoides_vcf_pruned_no_rel_0.4_pcadapt.bim"

## Create gds file
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "plink_pruned_no_rel_pcadapt.gds", cvt.chr="char")
snpgdsSummary("plink_pruned_no_rel_pcadapt.gds")
## Open the gds file
genofile <- snpgdsOpen("plink_pruned_no_rel_pcadapt.gds")
#The total number of samples: 13 
#The total number of SNPs: 5808 
###

## Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# metadata 
meta <- read.csv("Metadata.csv")
# sample names match
meta <- meta[match(sample.id, meta$id), ]

geno <- snpgdsGetGeno(genofile)  # matrix of SNPs x samples

#Convert 0/1/2 to hierfstat genotype format
geno.df <- as.data.frame(t(geno))

# Convert 0/1/2 to hierfstat genotype format
geno.df.hierf <- apply(geno.df, 2, function(x) {
  y <- rep(NA, length(x))
  y[x == 0] <- 11
  y[x == 1] <- 12
  y[x == 2] <- 22
  return(y)
})
geno.df.hierfstat <- as.data.frame(t(geno.df.hierf))

geno.df.hierfstat$pop <- meta$without.relatives
# hierfstat wants pop column to be first
geno.df.hierfstat <- geno.df.hierfstat[, c(ncol(geno.df.hierfstat), 1:(ncol(geno.df.hierfstat) - 1))]

genotype <- geno.df.hierfstat
genotype.AA <- apply(t(genotype), 2, function(x) {
  y <- rep(NA, length(x))
  y[x == 11] <- ("AA")
  y[x == 12] <- ("AB")
  y[x == 22] <- ("BB")
  return(y)
})
genotype <- as.data.frame(t(genotype.AA))

genotype$Sample <- as.factor(sample.id)
genotype$Ploidy <-2
genotype$Format <- "BiA"
genotype$Pop <- meta$without.relatives
genotype <- genotype[!is.na(genotype$Pop), ]

genotype <- genotype[,c("Sample", "Pop", "Ploidy", "Format", setdiff(names(genotype), c("Sample", "Pop", "Ploidy", "Format") ))]
cols <- 4:48477
genotype[] <-lapply(genotype, as.factor)
genotype$Ploidy <- 2
genotype.st <- stamppConvert(genotype, "r")

#### Fst ####

genotype.fst <- stamppFst(genotype.st,  nboots = 1000, percent = 95, nclusters=3)
genotype.fst[["Pvalues"]]
genotype.fst[["Fsts"]]

# confidence intervals
boot_data <- genotype.fst[["Bootstraps"]][, -c(1, 2)]
ci_lower <- apply(boot_data, 1, quantile, probs = 0.025, na.rm = TRUE)
ci_upper <- apply(boot_data, 1, quantile, probs = 0.975, na.rm = TRUE)

fst_summary <- data.frame(
  Pop1 = genotype.fst[["Bootstraps"]]$Population1,
  Pop2 = genotype.fst[["Bootstraps"]]$Population2,
  Fst = genotype.fst[["Fsts"]][lower.tri(genotype.fst[["Fsts"]])], # Extracts lower triangle values
  L95 = ci_lower,
  U95 = ci_upper,
  Pval = genotype.fst[["Pvalues"]][lower.tri(genotype.fst[["Pvalues"]])]
)

print(fst_summary)


#### Fst plot no rel pcadapt ####

# Prepare data
fst_mat  <- genotype.fst[["Fsts"]]
pval_mat <- genotype.fst[["Pvalues"]]
boot_data <- genotype.fst[["Bootstraps"]][, -c(1, 2)]

# Calculate CIs
ci_low <- apply(boot_data, 1, quantile, probs = 0.025, na.rm = TRUE)
ci_upp <- apply(boot_data, 1, quantile, probs = 0.975, na.rm = TRUE)

fst_df <- expand.grid(Pop1 = rownames(fst_mat), 
                      Pop2 = colnames(fst_mat), 
                      stringsAsFactors = FALSE)

# Add values and labels
fst_df$label   <- ""
fst_df$color_val <- NA_real_
boot_idx <- 1

for(i in 1:nrow(fst_df)) {
  p1 <- fst_df$Pop1[i]
  p2 <- fst_df$Pop2[i]
  
  r <- which(rownames(fst_mat) == p1)
  c <- which(colnames(fst_mat) == p2)
  
  if (r > c) {
    # LOWER TRIANGLE: Fst [CI]
    fst_val <- fst_mat[r, c]
    fst_df$label[i] <- sprintf("%.2f\n[%.2f, %.2f]", fst_val, ci_low[boot_idx], ci_upp[boot_idx])
    fst_df$color_val[i] <- fst_val
    boot_idx <- boot_idx + 1
  } else if (r < c) {
    # UPPER TRIANGLE: P-values 
    pv <- pval_mat[c, r] 
    fst_df$label[i] <- ifelse(is.na(pv), "", 
                              ifelse(pv < 0.001, "p < 0.001", sprintf("p = %.3f", pv)))
  }
}


# Final plot 
fst_no_rel_pcadapt <- ggplot(fst_df, aes(x = Pop2, y = Pop1)) +
  geom_tile(fill = "white", color = "grey90", lwd = 0.3) +
  geom_tile(aes(fill = color_val), color = "grey90", lwd = 0.3) +
  
  geom_text(
    aes(label = label),
    size = pt_to_mm(6),
    color = "black",
    lineheight = 0.75
  ) +
  
  scale_fill_gradientn(
    colors = c("#e0f3f8", "#abd9e9", "#74add1", "#4575b4"),
    name = expression(F[ST]),
    na.value = "transparent",
    guide = guide_colorbar(
      barwidth = unit(3, "mm"),
      barheight = unit(22, "mm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  
  scale_y_discrete(limits = rev(rownames(fst_mat))) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    
    axis.text = element_text(
      face = "bold",
      color = "black",
      size = 6
    ),
    
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.width = unit(3, "mm"),
    legend.key.height = unit(4, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, -4),
    
    plot.margin = margin(3, 1, 3, 3, unit = "pt")
  )
fst_no_rel_pcadapt

saveRDS(fst_no_rel_pcadapt, "fig_fst_no_rel_pcadapt.rds")
fst_no_rel_pcadapt <- readRDS("fig_fst_no_rel_pcadapt.rds")

#### PCA_no_rel pcadapt ####
genotype.st  -> genotype.st2
geno.nas <- genotype.st2[, !(names(genotype.st2) %in% c("Sample", "Pop", "pop.num", "ploidy", "format"))]
pop.vector <- as.factor(genotype.st2[[2]])
convert_to_genind_format <- function(x) {
  ifelse(is.na(x), NA,
         ifelse(x == 1, "11",
                ifelse(x == 0.5, "12",
                       ifelse(x == 0, "22", NA))))
}

geno.char <- as.data.frame(lapply(geno.nas, convert_to_genind_format), stringsAsFactors = FALSE)
geno.obj <- df2genind(geno.char, pop = pop.vector, ploidy = 2, NA.char = NA, sep = "")
dist.mtx <- dist(geno.obj)
pca <- dudi.pca(tab(geno.obj, NA.method = "mean"), scannf = F, nf = 2)
s.class(pca$li, fac = pop.vector, col = rainbow(length(unique(pop.vector))))


scores <- as.data.frame(pca$li)
scores$SampleID <- genotype.st2$Sample 
scores$Pop <- pop.vector

# Variance explained
eig <- 100 * pca$eig / sum(pca$eig)

# Calculate centroids per population
centroids <- scores %>%
  group_by(Pop) %>%
  summarise(
    Axis1 = mean(Axis1),
    Axis2 = mean(Axis2)
  )

# Add centroid coords to scores
scores <- scores %>%
  left_join(centroids, by = "Pop", suffix = c("", ".centroid"))

scores$label <- gsub("^\\./|\\-.*$", "", scores$SampleID)


# Plot
pca_plot_no_rel_pcadapt <- ggplot(scores, aes(x = Axis1, y = Axis2, color = Pop)) +
  # segments from sample to centroid
  geom_segment(aes(xend = Axis1.centroid, yend = Axis2.centroid),
               alpha = 0.4, linewidth = 0.3) +
  
  # centroid labels
  geom_text(data = centroids, aes(x = Axis1, y = Axis2, label = Pop, color = Pop),
            size = pt_to_mm(8), show.legend = FALSE, alpha = 0.6, fontface = "bold") +
  
  # points for samples
  geom_point(size = 1.5) +
  
  #  Sample numbers 
  geom_text_repel(aes(label = label), 
                  size = pt_to_mm(6), 
                  segment.size = 0.28,
                  segment.alpha = 0.6,
                  show.legend = FALSE,
                  max.overlaps = Inf) + 
  
  labs(
    x = paste0("PC1 (", round(eig[1], 1), "%)"),
    y = paste0("PC2 (", round(eig[2], 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "none")

pca_plot_no_rel_pcadapt


pca_plot_no_rel_pcadapt
saveRDS(pca_plot_no_rel_pcadapt, "fig_pca_no_rel_pcadapt.rds")
pca_plot_no_rel_pcadapt <- readRDS("fig_pca_no_rel_pcadapt.rds")


#### IBS dendrogram and tree no relatives pcadapt####
# Identity by State (IBS) Analysis
ibs <- snpgdsIBS(genofile, num.thread=2)
set.seed(100)
ibs.hc <- snpgdsHCluster(ibs)
rv <- snpgdsCutTree(ibs.hc)

dist_matrix <- as.dist(1 - ibs$ibs)
standard_hclust <- hclust(dist_matrix, method = "average")
standard_hclust$labels <- ibs$sample.id
ibs_phylo <- as.phylo(standard_hclust)
ibs_phylo$tip.label <- gsub("^\\./|\\-.*$", "", ibs_phylo$tip.label)
meta_subset <- meta[, c("id", "all")] 
colnames(meta_subset)[1] <- "label"
meta_subset$label <- gsub("^\\./|\\-.*$", "", meta_subset$label)

ibs_no_rel_pcadapt <- ggtree(ibs_phylo, layout="rectangular", linewidth = 0.3) %<+% meta_subset +
  geom_tiplab(size=pt_to_mm(6), offset=0.01) + 
  geom_tippoint(aes(color=all), size=1.5) +
  theme_tree2() +
  labs(color = "Group")+
  xlim(0, 0.19)

ibs_no_rel_pcadapt
saveRDS(ibs_no_rel_pcadapt, "fig_ibs_no_rel_pcadapt.rds")
ibs_no_rel_pcadapt <- readRDS("fig_ibs_no_rel_pcadapt.rds")

#### final plotting ####
kinship <- readRDS("fig_kinship.rds")
ibs_no_rel_pcadapt <- readRDS("fig_ibs_no_rel_pcadapt.rds")
pca_plot_no_rel_pcadapt <- readRDS("fig_pca_no_rel_pcadapt.rds")
fst_no_rel_pcadapt <- readRDS("fig_fst_no_rel_pcadapt.rds")
ibs_no_rel <- readRDS("fig_ibs_no_rel.rds")
pca_plot_no_rel <- readRDS("fig_pca_no_rel.rds")
fst_no_rel <- readRDS("fig_fst_no_rel.rds")
ibs_no_art <- readRDS("fig_ibs_no_art.rds")
pca_plot_no_art <- readRDS("fig_pca_no_art.rds")

#### supplementary fig ####
fig_theme_sup <- theme(
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6),
  legend.title = element_text(size = 6),
  legend.text = element_text(size = 6),
  strip.text = element_text(size = 6),
  text = element_text(size = 6)
)


pca_colors <- c(
  "D.CC" = "#00C19A", 
  "D.MF" = "#00A9FF", 
  "S.CC" = "#ED68ED", 
  "S.MF" = "#FF68A1"
)
pca_plot_no_art2 <- pca_plot_no_art +
  scale_color_manual(values = pca_colors) +
  fig_theme_sup

pca_plot_no_rel2 <- pca_plot_no_rel +
  scale_color_manual(values = pca_colors) +
  fig_theme_sup

pca_plot_no_rel_pcadapt2 <- pca_plot_no_rel_pcadapt +
  scale_color_manual(values = pca_colors) +
  fig_theme_sup

ibs_no_rel2 <- ibs_no_rel +
  scale_color_manual(values = pca_colors) +
  fig_theme_sup

ibs_no_rel_pcadapt2 <- ibs_no_rel_pcadapt +
  scale_color_manual(values = pca_colors) +
  fig_theme_sup

fst_no_rel_pcadapt2 <- fst_no_rel_pcadapt+
  fig_theme_sup 

combined_sup <- (
  (pca_plot_no_art2 | pca_plot_no_rel2 | pca_plot_no_rel_pcadapt2) /
    (fst_no_rel_pcadapt2 | ibs_no_rel2 | ibs_no_rel_pcadapt2)
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.tag = element_text(size = 9, face = "bold")
    ))


ggsave(
  "/home/gospozha/haifa/cayman/manuscript/revision/figs/combined_figS3.pdf",
  plot = combined_sup,
  width = 18.4,
  height = 10,
  units = "cm",
  device = "pdf",
  useDingbats = FALSE
)

#### fig 2 ####
ibs_no_art_tagged <- ibs_no_art +
  ggplot2::labs(tag = "A") +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 9)
  )+ scale_color_manual(values = pca_colors) +
  fig_theme_sup
  

kinship_tagged <- kinship +
  ggplot2::labs(tag = "B") +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 9)
  )+ scale_color_manual(values = pca_colors) 

fst_no_rel_tagged <- fst_no_rel +
  ggplot2::labs(tag = "C") +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 9)
  )+ scale_color_manual(values = pca_colors) +
  fig_theme_sup

combined_fig2 <- (ibs_no_art_tagged | (kinship_tagged / fst_no_rel_tagged)) +
  patchwork::plot_layout(widths = c(1.3, 0.7))

ggsave(
  "/home/gospozha/haifa/cayman/manuscript/revision/figs/combined_fig2.pdf",
  plot = combined_fig2,
  width = 18.4,
  height = 10,
  units = "cm",
  device = "pdf",
  useDingbats = FALSE
)

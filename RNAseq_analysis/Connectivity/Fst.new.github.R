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

setwd("~/haifa/cayman/rna/connect/admixture/new/")

#### no artificial clones ####
# all except for the artificial clones and samples below missing 0.9 threshold
## Load the pruned SNPs

bed.fn <- "P_astreoides_vcf_pruned_rel_all.bed"
fam.fn <- "P_astreoides_vcf_pruned_rel_all.fam"
bim.fn <- "P_astreoides_vcf_pruned_rel_all.bim"

## Create gds file
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "plink_pruned_rel_all.gds", cvt.chr="char")
snpgdsSummary("plink_pruned_rel_all.gds")
## Open the gds file
genofile <- snpgdsOpen("plink_pruned_rel_all.gds")

## Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Load your metadata (assumed to be a CSV with columns: sample, site, depth)
meta <- read.csv("../../Metadata.csv")
# Ensure sample names match
meta <- meta[match(sample.id, meta$id), ]

geno <- snpgdsGetGeno(genofile)  # matrix of SNPs x samples

# hierfstat expects genotypes in 1/2 digit allele coding (e.g., 11, 12, 22), not 0/1/2
# We can fake this by converting:

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

# convert to StAMPP format
genotype <- geno.df.hierfstat
genotype.AA <- apply(t(genotype), 2, function(x) {
  y <- rep(NA, length(x))
  y[x == 11] <- ("AA")
  y[x == 12] <- ("AB")
  y[x == 22] <- ("BB")
  return(y)
})
genotype <- as.data.frame(t(genotype.AA))

#genotype$Pop <- as.factor(geno.df.hierfstat$pop)
genotype$Sample <- as.factor(sample.id)
genotype$Ploidy <-2
genotype$Format <- "BiA"
genotype$Pop <- meta$without.relatives
#genotype$Pop <- c("D", "D", "S", "D", "S", "D", "S", "D", "D", "S", "D", "D", "S")
genotype <- genotype[!is.na(genotype$Pop), ]

genotype <- genotype[,c("Sample", "Pop", "Ploidy", "Format", setdiff(names(genotype), c("Sample", "Pop", "Ploidy", "Format") ))]
cols <- 4:48477
genotype[] <-lapply(genotype, as.factor)
genotype$Ploidy <- 2
genotype.st <- stamppConvert(genotype, "r")

genotype.fst <- stamppFst(genotype.st,  nboots = 1000, percent = 95, nclusters=3)
genotype.fst[["Pvalues"]]
genotype.fst[["Fsts"]]

#### plot ####
fst <- genotype.fst[["Fsts"]]
# Apply mask to upper triangle
mask_upper_triangle <- function(mat) {
  mat[upper.tri(mat)] <- NA
  return(mat)
}
fst_lower <- mask_upper_triangle(fst)

# Create label matrix
#label_matrix <- matrix(sprintf("%.2f", fst), nrow = 3)
label_matrix <- matrix(sprintf("%.2f", fst), nrow = 4)
label_matrix[upper.tri(label_matrix)] <- NA
label_matrix[is.na(label_matrix)] <- ""

# Clustering based on symmetric version (NAs replaced with 0 for distance calc)
fst_clust <- fst
fst_clust[is.na(fst_clust)] <- 0
row_order <- hclust(dist(fst_clust), method = "complete")

# Get dimensions
n_rows <- nrow(fst_lower)
n_cols <- ncol(fst_lower)

# Create 5 breaks and 5 matching colors
breaks <- pretty(range(fst, na.rm = TRUE), n = 5)
colors <- colorRampPalette(c("#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695"))(length(breaks))

col_fun <- colorRamp2(breaks = breaks, colors = colors)
custom_order <- c("MF.S", "MF.D", "CC.D", "CC.S")
#custom_order <- c("MF.D", "CC.D", "CC.S")
# Plot
ht <- Heatmap(
  fst_lower,
  name = "Fst",
  col = col_fun,
  na_col = "white",
  cluster_rows = FALSE,
  #column_order = custom_order,
  cluster_columns = FALSE,
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  rect_gp = gpar(col = NA),
  cell_fun = function(j, i, x, y, width, height, fill) {
    label <- label_matrix[i, j]
    
    # Draw bottom border for last row
    if (i == n_rows) {
      grid.lines(x = unit.c(x - 0.5 * width, x + 0.5 * width),
                 y = unit.c(y - 0.5 * height, y - 0.5 * height),
                 gp = gpar(col = "black", lwd = 1))
    }
    # Draw left border for first column
    if (j == 1) {
      grid.lines(x = unit.c(x - 0.5 * width, x - 0.5 * width),
                 y = unit.c(y - 0.5 * height, y + 0.5 * height),
                 gp = gpar(col = "black", lwd = 1))
    }
    
    # Add text label
    if (label != "") {
      grid.text(label, x, y, gp = gpar(fontsize = 9, col = "black"))
    }
  },
  heatmap_legend_param = list(
    at = pretty(range(fst, na.rm = TRUE), n = 5),
    labels = as.character(pretty(range(fst, na.rm = TRUE), n = 5)),
    title = "Fst",
    title_gp = gpar(fontsize = 11, fontface = "plain"),
    labels_gp = gpar(fontsize = 9)
  )
)
ht
draw(ht)

#### join with kinship ####

# A) Kinship network saved as ggplot
p1 <- readRDS("~/haifa/cayman/rna/connect/admixture/fig_kinship.rds")

# If you want plot margins, do it on the ggplot BEFORE converting to grob:
#p1 <- p1 + theme(plot.margin = margin(t = 10, r = 5, b = 5, l = 5))
g1 <- ggplotGrob(p1)

# B) ComplexHeatmap -> capture as grob
g2 <- grid.grabExpr(
  draw(ht, heatmap_legend_side = "right")
)

# Add panel labels (move a bit right with x = 0.06 npc)
g1_t <- arrangeGrob(
  g1,
  top = textGrob("A",
                 x = unit(0.02, "npc"), just = "left",
                 gp = gpar(fontsize = 12))
)

g2_t <- arrangeGrob(
  g2,
  top = textGrob("B",
                 x = unit(0.02, "npc"), just = "left",
                 gp = gpar(fontsize = 12))
)

# Combine with a wider gap between panels
combined <- arrangeGrob(
  g1_t,
  nullGrob(),   # spacer column
  g2_t,
  ncol = 3,
  widths = unit.c(
    unit(1, "null"),
    unit(12, "mm"),  # GAP between A and B (increase to 15–20 mm if needed)
    unit(1, "null")
  )
)

ggsave("Fst.pdf", combined, width = 200, height = 110, units = "mm")


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


# Extract individual scores (coordinates)
scores <- as.data.frame(pca$li)
scores$Pop <- pop.vector
scores$SampleID <- rownames(scores)

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

# Plot
pca_plot2 <- ggplot(scores, aes(x = Axis1, y = Axis2, color = Pop)) +
  # segments from sample to centroid
  geom_segment(aes(xend = Axis1.centroid, yend = Axis2.centroid),
               alpha = 0.4, linewidth = 0.3) +
  # centroid labels (drawn first, so underneath points)
  geom_text(data = centroids, aes(x = Axis1, y = Axis2, label = Pop, color = Pop),
            , size = 6, show.legend = FALSE, alpha = 0.6) +
  # points for samples (drawn after, so on top)
  geom_point(size = 3) +
  # centroid markers (optional cross shape)
  #geom_point(data = centroids, aes(x = Axis1, y = Axis2, color = Pop),
  #           size = 1, shape = 4, stroke = 1.5) +
  labs(
    x = paste0("PC1 (", round(eig[1], 1), "%)"),
    y = paste0("PC2 (", round(eig[2], 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "none")

saveRDS(pca_plot2, "pca_plot_all.RDS")

percent = pca$eig/sum(pca$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,12),
        names.arg = round(percent, 1))


png("sclass.png", width=1200, height=1200, res=150)
s.class(pca$li, fac = pop.vector, col = rainbow(length(unique(pop.vector))))
dev.off()


# plotting

ht_grob <- grid.grabExpr(draw(ht))

plotA <- ggdraw() +
  draw_label("A. Heatmap of pairwise Fst", fontface = "bold", size = 12,
             x = 0.5, y = 0.98, hjust = 0.5) +
  draw_plot(ht_grob, 0, 0, 1, 0.93)


plotB <- ggdraw() +
  draw_label("B. PCA plot", fontface = "bold", size = 12,
             x = 0.5, y = 0.98, hjust = 0.5) +
  draw_plot(combined, 0, 0, 1, 0.93)

final <- plot_grid(plotA, plotB, ncol = 1, rel_widths = c(1, 1))

final
ggsave("combined_all_clones.jpg", final, width = 10, height = 5)


#### heterozygosy analysis ####

# Make sure all columns are numeric (except pop)
geno.df.hierfstat[,-1] <- lapply(geno.df.hierfstat[,-1], function(x) as.numeric(as.character(x)))
geno.df.hierfstat$pop <- as.factor(geno.df.hierfstat$pop)

# Convert to hierfstat format
gen_hier <- geno.df.hierfstat

# Check structure
str(gen_hier)
# First column = pop (factor), remaining columns = numeric genotypes 11/12/22

#  basic statistics per population 
bs <- basic.stats(gen_hier)

# Observed heterozygosity Ho per population
Ho <- colMeans(bs$Ho, na.rm=TRUE)

# Expected heterozygosity He per population
He <- colMeans(bs$Hs, na.rm=TRUE)

# FIS per population
Fis <- colMeans(bs$Fis, na.rm=TRUE)

# allelic richness per population 
ar <- allelic.richness(gen_hier)


# overall allelic richness per population - the mean across loci
allelic_richness_pop <- colMeans(ar$Ar, na.rm=TRUE)

results <- data.frame(
  Site = names(allelic_richness_pop),  # population names
  Ho = Ho,
  He = He,
  FIS = Fis,
  AllelicRichness = allelic_richness_pop
)


print(results)
write.csv(results, "stat.rel.all.csv")


#### first degree only in MF.S ####

## Load the pruned SNPs

bed.fn <- "P_astreoides_vcf_pruned_rel_mfs.bed"
fam.fn <- "P_astreoides_vcf_pruned_rel_mfs.fam"
bim.fn <- "P_astreoides_vcf_pruned_rel_mfs.bim"

## Create gds file
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "plink_pruned_rel_mfs2.gds", cvt.chr="char")
snpgdsSummary("plink_pruned_rel_mfs2.gds")
## Open the gds file
genofile <- snpgdsOpen("plink_pruned_rel_mfs2.gds")

## Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Load your metadata (assumed to be a CSV with columns: sample, site, depth)
meta <- read.csv("../../Metadata.csv")
# Ensure sample names match
meta <- meta[match(sample.id, meta$id), ]

geno <- snpgdsGetGeno(genofile)  # matrix of SNPs x samples

# hierfstat expects genotypes in 1/2 digit allele coding (e.g., 11, 12, 22), not 0/1/2
# We can fake this by converting:

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

# convert to StAMPP format
genotype <- geno.df.hierfstat
genotype.AA <- apply(t(genotype), 2, function(x) {
  y <- rep(NA, length(x))
  y[x == 11] <- ("AA")
  y[x == 12] <- ("AB")
  y[x == 22] <- ("BB")
  return(y)
})
genotype <- as.data.frame(t(genotype.AA))

#genotype$Pop <- as.factor(geno.df.hierfstat$pop)
genotype$Sample <- as.factor(sample.id)
genotype$Ploidy <-2
genotype$Format <- "BiA"
genotype$Pop <- meta$without.relatives
#genotype$Pop <- c("D", "D", "S", "D", "S", "D", "S", "D", "D", "S", "D", "D", "S")
genotype <- genotype[!is.na(genotype$Pop), ]

genotype <- genotype[,c("Sample", "Pop", "Ploidy", "Format", setdiff(names(genotype), c("Sample", "Pop", "Ploidy", "Format") ))]
cols <- 4:48477
genotype[] <-lapply(genotype, as.factor)
genotype$Ploidy <- 2
genotype.st <- stamppConvert(genotype, "r")

genotype.fst <- stamppFst(genotype.st,  nboots = 1000, percent = 95, nclusters=3)
genotype.fst[["Pvalues"]]
genotype.fst[["Fsts"]]

#### plot ####
fst <- genotype.fst[["Fsts"]]
# Apply mask to upper triangle
mask_upper_triangle <- function(mat) {
  mat[upper.tri(mat)] <- NA
  return(mat)
}
fst_lower <- mask_upper_triangle(fst)

# Create label matrix
#label_matrix <- matrix(sprintf("%.2f", fst), nrow = 3)
label_matrix <- matrix(sprintf("%.2f", fst), nrow = 4)
label_matrix[upper.tri(label_matrix)] <- NA
label_matrix[is.na(label_matrix)] <- ""

# Clustering based on symmetric version (NAs replaced with 0 for distance calc)
fst_clust <- fst
fst_clust[is.na(fst_clust)] <- 0
row_order <- hclust(dist(fst_clust), method = "complete")

# Get dimensions
n_rows <- nrow(fst_lower)
n_cols <- ncol(fst_lower)

# Create 5 breaks and 5 matching colors
breaks <- pretty(range(fst, na.rm = TRUE), n = 5)
colors <- colorRampPalette(c("#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695"))(length(breaks))

col_fun <- colorRamp2(breaks = breaks, colors = colors)
custom_order <- c("MF.S", "MF.D", "CC.D", "CC.S")
#custom_order <- c("MF.D", "CC.D", "CC.S")
# Plot
ht <- Heatmap(
  fst_lower,
  name = "Fst",
  col = col_fun,
  na_col = "white",
  cluster_rows = FALSE,
  #column_order = custom_order,
  cluster_columns = FALSE,
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  rect_gp = gpar(col = NA),
  cell_fun = function(j, i, x, y, width, height, fill) {
    label <- label_matrix[i, j]
    
    # Draw bottom border for last row
    if (i == n_rows) {
      grid.lines(x = unit.c(x - 0.5 * width, x + 0.5 * width),
                 y = unit.c(y - 0.5 * height, y - 0.5 * height),
                 gp = gpar(col = "black", lwd = 1))
    }
    # Draw left border for first column
    if (j == 1) {
      grid.lines(x = unit.c(x - 0.5 * width, x - 0.5 * width),
                 y = unit.c(y - 0.5 * height, y + 0.5 * height),
                 gp = gpar(col = "black", lwd = 1))
    }
    
    # Add text label
    if (label != "") {
      grid.text(label, x, y, gp = gpar(fontsize = 9, col = "black"))
    }
  },
  heatmap_legend_param = list(
    at = pretty(range(fst, na.rm = TRUE), n = 5),
    labels = as.character(pretty(range(fst, na.rm = TRUE), n = 5)),
    title = "Fst",
    title_gp = gpar(fontsize = 11, fontface = "plain"),
    labels_gp = gpar(fontsize = 9)
  )
)
ht
draw(ht)

#### heterozygousy analysis ####

# Make sure all columns are numeric (except pop)
geno.df.hierfstat[,-1] <- lapply(geno.df.hierfstat[,-1], function(x) as.numeric(as.character(x)))
geno.df.hierfstat$pop <- as.factor(geno.df.hierfstat$pop)

# Convert to hierfstat format
gen_hier <- geno.df.hierfstat

# Check structure
str(gen_hier)
# First column = pop (factor), remaining columns = numeric genotypes 11/12/22

#  basic statistics per population 
bs <- basic.stats(gen_hier)

# Observed heterozygosity Ho per population
Ho <- colMeans(bs$Ho, na.rm=TRUE)

# Expected heterozygosity He per population
He <- colMeans(bs$Hs, na.rm=TRUE)

# FIS per population
Fis <- colMeans(bs$Fis, na.rm=TRUE)

# allelic richness per population
ar <- allelic.richness(gen_hier)

# Check structure of allelic richness object
str(ar$Ar)
# Typically it's: loci x populations

# overall allelic richness per population - the mean across loci
allelic_richness_pop <- colMeans(ar$Ar, na.rm=TRUE)

# Now combine into results table
results <- data.frame(
  Site = names(allelic_richness_pop),  # population names
  Ho = Ho,
  He = He,
  FIS = Fis,
  AllelicRichness = allelic_richness_pop
)


print(results)
write.csv(results, "stat.rel.mfsonly.csv")

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


# Extract individual scores (coordinates)
scores <- as.data.frame(pca$li)
scores$Pop <- pop.vector
scores$SampleID <- rownames(scores)

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

# Plot
pca_plot3 <- ggplot(scores, aes(x = Axis1, y = Axis2, color = Pop)) +
  # segments from sample to centroid
  geom_segment(aes(xend = Axis1.centroid, yend = Axis2.centroid),
               alpha = 0.4, linewidth = 0.3) +
  # centroid labels (drawn first, so underneath points)
  geom_text(data = centroids, aes(x = Axis1, y = Axis2, label = Pop, color = Pop),
            , size = 6, show.legend = FALSE, alpha = 0.6) +
  # points for samples (drawn after, so on top)
  geom_point(size = 3) +
  # centroid markers (optional cross shape)
  #geom_point(data = centroids, aes(x = Axis1, y = Axis2, color = Pop),
  #           size = 1, shape = 4, stroke = 1.5) +
  labs(
    x = paste0("PC1 (", round(eig[1], 1), "%)"),
    y = paste0("PC2 (", round(eig[2], 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "none")

saveRDS(pca_plot3, "pca_plot_MFonly.RDS")


#### combine PCA plots ####

combined <- (pca_plot2 | pca_plot3) + 
  plot_annotation(tag_levels = 'A')

combined
ggsave("~/haifa/cayman/rna/mapping/github/fin/combined_PCA.jpg", combined, width = 8, height = 4)
ggsave("~/haifa/cayman/rna/mapping/github/fin/combined_PCA.pdf", combined, width = 8, height = 4)

# with a new PCA plot

ht_grob <- grid.grabExpr(draw(ht))

plotA <- ggdraw() +
  draw_plot(ht_grob, 0, 0, 1, 0.93)


plotB <- ggdraw() +
  draw_plot(combined, 0, 0, 1, 0.93)

final <- plot_grid(plotA, plotB, ncol = 1, rel_widths = c(1, 1))

final
ggsave("~/haifa/cayman/rna/mapping/github/fin/combined_PCA_Fst.jpg", final, width = 6, height = 6)
ggsave("~/haifa/cayman/rna/mapping/github/fin/combined_PCA_Fst.pdf", final, width = 6, height = 6)



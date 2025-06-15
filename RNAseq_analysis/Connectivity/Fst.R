library(gdsfmt)
library(SNPRelate)
library(hierfstat)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)


#setwd("/lustre1/home/mass/eskalon/Porites/analysis/snp/Fst")
#setwd("~/haifa/cayman/rna/connect/")

load("./server/Fst3.RData")

## Load the pruned SNPs

bed.fn <- "P_astreoides_vcf_pruned.bed"
fam.fn <- "P_astreoides_vcf_pruned.fam"
bim.fn <- "P_astreoides_vcf_pruned.bim"

## Create gds file
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "plink_pruned.gds")
snpgdsSummary("plink_pruned.gds")

# The file name: /home/gospozha/haifa/cayman/rna/connect/plink_pruned.gds 
# The total number of samples: 31 
# The total number of SNPs: 48473 
# SNP genotypes are stored in SNP-major mode (Sample X SNP).
# The number of valid samples: 31 

# The number of biallelic unique SNPs: 62 
## Open the gds file
genofile <- snpgdsOpen("plink_pruned.gds")

## Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Load your metadata (assumed to be a CSV with columns: sample, site, depth)
meta <- read.csv("Metadata.csv")
# Ensure sample names match
meta <- meta[match(sample.id, meta$id), ]

meta$group <- interaction(meta$site, meta$depth, drop = TRUE)
group.factor <- meta$group

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
geno.df.hierfstat$pop <- group.factor
# hierfstat wants pop column to be first
geno.df.hierfstat <- geno.df.hierfstat[, c(ncol(geno.df.hierfstat), 1:(ncol(geno.df.hierfstat) - 1))]


keep_ids <- c("MF.S", "MF.D", "CC.S", "CC.D")

geno.df.hierfstat.cut <- geno.df.hierfstat %>% dplyr::filter(pop %in% keep_ids)


# For pairwise Fst
pairwise.fst <- pairwise.WCfst(geno.df.hierfstat.cut)
pairwise.fst
pairwise.fst[pairwise.fst < 0] <- 0
pairwise.fst
# Bootstrapping to test significance
boot.fst <- boot.ppfst(dat = geno.df.hierfstat, nboot = 1000)
boot.fst
boot.fst$ll[boot.fst$ll < 0] <- 0
boot.fst$ul[boot.fst$ul < 0] <- 0


# Convert each matrix to long format (vector + names)
get_long_df <- function(mat, value_name) {
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("pop1", "pop2", value_name)
  return(df)
}

fst_df <- get_long_df(pairwise.fst, "fst")
ll_df  <- get_long_df(boot.fst$ll, "ci.low")
ul_df  <- get_long_df(boot.fst$ul, "ci.high")

# Merge all into one
df <- Reduce(function(x, y) merge(x, y, by = c("pop1", "pop2")), list(fst_df, ll_df, ul_df))

# Optional: Remove duplicates (keep only upper triangle)
df <- na.omit(df)

# Calculate CI width
df$ci.width = df$ci.high - df$ci.low

# Sort by Fst or CI width or both
df <- df[order(df$fst, decreasing = TRUE), ]
print(df)


save.image(file="Fst3.RData")

#### Plots ####

# Error bar plot
error <- ggplot(df, aes(x = reorder(paste(pop1, pop2, sep = " vs "), fst), y = fst)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci.low, ymax = ci.high), width = 0.1) +
  coord_flip() +
  theme_minimal() +
  ylab("Fst") + xlab("Population Comparison") +
  ggtitle("Pairwise Fst with 95% confidence intervals")
error
ggsave("error.fst3.jpg", error, width = 7, height = 7)


# Get all unique population names
all_pops <- sort(unique(c(df$pop1, df$pop2)))

# Create a full matrix of NA
label_matrix <- matrix(NA, nrow = length(all_pops), ncol = length(all_pops),
                       dimnames = list(all_pops, all_pops))

# Fill in the label matrix from fst_df
for (i in seq_len(nrow(df))) {
  pop1 <- df$pop1[i]
  pop2 <- df$pop2[i]
  label <- sprintf("%.2f\n[%.2f–%.2f]", df$fst[i], df$ci.low[i], df$ci.high[i])
  label_matrix[pop1, pop2] <- label
  label_matrix[pop2, pop1] <- label  # optional: mirror for symmetry
}

# Reorder both matrices to match
ord <- rownames(pairwise.fst)  # assumes this is already ordered well
label_matrix <- label_matrix[ord, ord]

# Mask upper triangle
mask_upper_triangle <- function(mat) {
  mat[upper.tri(mat)] <- NA
  return(mat)
}

# Apply mask
fst_lower <- mask_upper_triangle(pairwise.fst)
labels_lower <- mask_upper_triangle(label_matrix)
labels_lower[is.na(labels_lower)] <- ""

# Color scale
col_fun <- colorRamp2(
  breaks = pretty(range(pairwise.fst, na.rm = TRUE), n = 5),
  colors = c("#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
)

# Dendrogram for rows only
row_order <- hclust(dist(pairwise.fst), method = "complete")

# Get number of rows and columns for framing
n_rows <- nrow(fst_lower)
n_cols <- ncol(fst_lower)

# Heatmap
jpeg("fst_heatmap.jpg", width = 2000, height = 2000, res = 300)
ht <- Heatmap(
  fst_lower,
  name = "Fst",
  col = col_fun,
  na_col = "white",
  cluster_rows = as.dendrogram(row_order),
  column_names_rot = 0,
  cluster_columns = FALSE,
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_side = "bottom",
  rect_gp = gpar(col = NA),  # No default grid lines
  cell_fun = function(j, i, x, y, width, height, fill) {
    label <- labels_lower[i, j]
    
    # Always draw bottom and left border for the outermost row/column
    if (i == n_rows) {
      grid.lines(x = unit.c(x - 0.5 * width, x + 0.5 * width),
                 y = unit.c(y - 0.5 * height, y - 0.5 * height),
                 gp = gpar(col = "black", lwd = 1))
    }
    if (j == 1) {
      grid.lines(x = unit.c(x - 0.5 * width, x - 0.5 * width),
                 y = unit.c(y - 0.5 * height, y + 0.5 * height),
                 gp = gpar(col = "black", lwd = 1))
    }
    
    # Add label only for non-NA cells
    if (!is.na(fst_lower[i, j]) && label != "") {
      grid.text(label, x, y, gp = gpar(fontsize = 12, col = "black"))
    }
  },
  heatmap_legend_param = list(
    title = expression(F[ST]),
    at = pretty(range(pairwise.fst, na.rm = TRUE), n = 5),
    labels = as.character(pretty(range(pairwise.fst, na.rm = TRUE), n = 5)),
    title_gp = gpar(fontsize = 16, fontface = "bold"),  # Title size and style
    labels_gp = gpar(fontsize = 12) 
  )
)

draw(ht, column_title = "Pairwise Fst with 95% confidence intervals", column_title_gp = gpar(fontsize = 17))
dev.off()

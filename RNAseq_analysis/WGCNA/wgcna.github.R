library("tidyverse")
library("DESeq2")
library("RColorBrewer")
library("WGCNA")
library("flashClust")
library("gridExtra")
library("ComplexHeatmap")
library("goseq")
library("dplyr")
library("clusterProfiler")
library("simplifyEnrichment")
library(dendsort)
library(edgeR)
library(multcompView)
library(ggforce)
library(gridExtra)
library(tibble)
library(ggpattern)


load("/home/gospozha/haifa/cayman/rna/mapping/github/rin/070725.vst.RData")

# setting working directory 
setwd("/home/gospozha/haifa/cayman/rna/mapping/github")

#### reading files ####
countData  <- read.csv2('CountMatrix.csv', header=TRUE, row.names=1, sep=',', check.names = FALSE)
MetaData <- read.csv2('Metadata.origin.csv', header=TRUE, sep=',')

MetaData$condition <- factor(MetaData$condition, levels = c("S", "D", "SS", "DD", "SD", "DS"))
MetaData$site <- factor(MetaData$site, levels = c("MF","CC"))
MetaData$translocation <- as.factor(MetaData$translocation)
MetaData$month <- as.factor(MetaData$month)
MetaData$origin <- as.factor(MetaData$origin)

# filtering genes 
smallestGroupSize <- 4
keep <- rowSums(countData >= 10) >= smallestGroupSize
countData <- countData[keep,]

#### vst transform (already done during DESeq2 analysis) ####
# vst transformation using deseq2
# creating DESeq2 object 
# dds <- DESeqDataSetFromMatrix(countData = countData,
#                               colData = MetaData,
#                               design = ~ 1)
# SF <- estimateSizeFactors(dds) 
# print(sizeFactors(SF))
# # variance-stabilizing transformation
# vst <- vst(dds)
# # since vst does not remove variation that can be associated with covariates, 
# # we manually remove the effect of covariates to be able to visualize it on PCA
# mat <- assay(vst)
# mm <- model.matrix(~condition, colData(vst))
# mat <- limma::removeBatchEffect(mat, batch=vst$site, design=mm)
# assay(vst) <- mat
#vsd <- assay(vst)

vsd <- read.csv(file="rin/GO/vsd.counts.csv", check.names = F, header = T, row.names = 1)

#### check for outliers ####
sampleTree <- hclust(dist(t(vsd)), method = "average")
# Plot the sample tree
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree,
     main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2
)
#Plot a line showing the cut-off
abline(h = 220, col = "red")

# remove two outliers since WGCNA is more sensitive to them
datExpr <- as.data.frame(t(vsd))
keepSamples <- !(rownames(datExpr) %in% c("78-MF-DS-54", "91-CC-DS-89"))  # Replace with actual IDs
datExpr <- datExpr[keepSamples, ]
MetaData <- MetaData[keepSamples, ]

#Check for genes and samples with too many missing values with goodSamplesGenes. There shouldn't be any because we performed pre-filtering
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK #Should return TRUE

#### pick soft power threshold ####
powers = c(c(1:10), seq(from = 9, to=25, by=1))
sft <- pickSoftThreshold(datExpr,
                        powerVector = powers,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)

sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()
# manual calculations
softPower=11 #w/o outliers

#### get modules ####
# automatic modules
bwnet <- blockwiseModules(datExpr,
                          maxBlockSize = 16000, # What size chunks (how many genes) the calculations should be run in
                          TOMType = "signed", # topological overlap matrix
                          power = softPower, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          randomSeed = 1234, # there's some randomness associated with this calculation
                          # so we should set a seed
                          minModuleSize = 30,
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          mergeCutHeight = 0.2,
                          deepSplit = 2
                          )

#readr::write_rds(bwnet,
                 file = file.path("wgcna_results.rin.rrna.0.2.2.vst.outliers.RDS"))
#bwnet <- read_rds("wgcna_results.rin.rrna.0.2.2.vst.outliers.RDS")

# n of genes in each module
table(bwnet$colors)
# 10000 genes are not assigned to any module (gray module)

# plotting dendrogram
mergedColors = labels2colors(bwnet$colors)
plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)

# We can also see if our eigengenes relate to our metadata labels. 
# First we double check that our samples are still in order.

all.equal(MetaData$id, rownames(module_eigengenes))

#### limma on modules ####
# Create the design matrix from the `time_point` variable
des_mat <- model.matrix(~ MetaData$condition)

# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

head(stats_df)

# boxplots of certain modules
module_26_df <- module_eigengenes %>%
  tibble::rownames_to_column("id") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(MetaData %>%
                      dplyr::select(id, condition),
                    by = c("id" = "id"))

ggplot(module_26_df,
  aes(
    x = condition,
    y = ME2,
    color = condition)) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()

# which genes make up modules
gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))

gene_module_key %>%
  dplyr::filter(module == "ME2")

##### function for a heatmap #####
make_module_heatmap <- function(module_name,
                                expression_mat = datExpr,
                                metadata_df = MetaData,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with refinebio_accession_code and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its refinebio_accession_code
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("id")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(id, condition) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "id") %>%
    # Arrange by patient and time point
    dplyr::arrange(condition) %>%
    # Store sample
    tibble::column_to_rownames("id")
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    condition = col_annot_df$condition,
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in time_point
    col = list(condition = c("S" = "#f1a340", "D" = "#998ec3", 
                             "SS" = "gray", "DD" = "violet", 
                             "SD" = "pink", "DS" = "lightblue"))
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = FALSE,
                                     show_column_names = FALSE
  )
  
  # Return heatmap
  return(heatmap)
}

mod_2_heatmap <- make_module_heatmap(module_name = "ME18")
mod_49_heatmap <- make_module_heatmap(module_name = "ME4")
mod_49_heatmap

#Relating modules to characteristics and identifying important genes
#Defining the number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# making binary traits vector
MetaData_bin <- MetaData %>%
  dplyr::select(id, condition)  %>%
  mutate(value = 1) %>%  # Create a helper column for filling binary values
  pivot_wider(names_from = condition, values_from = value, values_fill = list(value = 0)) %>%
  column_to_rownames("id")

#Recalculating MEs with label colors
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, MetaData_bin, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#sizeGrWindow(15,8)

#Displaying correlations and its p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
#par(mar = c(6, 8.5, 3, 3))

#Displaying the correlation values in a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MetaData_bin),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# Plot as clustered Heatmap
#add bold sigignificant p-values, dendrogram with WGCNA MEtree cut-off, module clusters

#Create list of pvalues for eigengene correlation with specific life stages
heatmappval <- signif(moduleTraitPvalue, 1)

#Make list of heatmap row colors
htmap.colors <- names(MEs)
htmap.colors <- gsub("ME", "", htmap.colors)

row_dend = dendsort(hclust(dist(moduleTraitCor)))
col_dend = dendsort(hclust(dist(t(moduleTraitCor))))

# clustering modules
#### complex heatmap ####
MEDiss = 1-cor(MEs)
METree = flashClust(as.dist(MEDiss), method = "average")
MEtreePlot = plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

ht=Heatmap(moduleTraitCor, name = "Eigengene", column_title = "Module-trait eigengene correlation", 
           col = blueWhiteRed(50), 
           row_names_side = "left", row_dend_side = "left",
           width = unit(3, "in"), height = unit(7, "in"), 
           column_order = 1:6, column_dend_reorder = FALSE, cluster_columns = hclust(dist(t(moduleTraitCor)), method = "average"), column_split = 6, column_dend_height = unit(0.5, "in"),
           cluster_rows = METree, row_split = 9, row_gap = unit(2.5, "mm"), border = TRUE,
           cell_fun = function(j, i, x, y, w, h, col) {
             if(heatmappval[i, j] <= 0.05) {
               grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
             }
             else {
               grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "plain"))
             }},
           column_names_gp =  gpar(fontsize = 10),
           row_names_gp = gpar(fontsize = 10, alpha = 0.75, border = TRUE, fill = htmap.colors))
ht <- draw(ht)



rcl.list <- row_order(ht)  #Extract clusters (output is a list)

clusters <- lapply(rcl.list, function(idx) rownames(moduleTraitCor)[idx])
names(clusters) <- seq_along(clusters)
moduleCluster <- stack(clusters)
colnames(moduleCluster) <- c("moduleColor", "moduleCluster")


#### performing DE analysis on modules clusters ####

#View module eigengene data and make dataframe for Strader plots.

head(MEs)
names(MEs)
Strader_MEs <- MEs
Strader_MEs$condition <- MetaData$condition
Strader_MEs$sample_id <- rownames(Strader_MEs)
head(Strader_MEs)

# Initialize output dataframe
cluster_list <- list()

for (i in sort(unique(moduleCluster$moduleCluster))) {
  colors <- moduleCluster %>%
    filter(moduleCluster == i) %>%
    pull(moduleColor)
  
  me_names <- paste0(colors)
  selected <- select(Strader_MEs, all_of(me_names))
  
  cluster_mean <- rowMeans(selected, na.rm = TRUE)
  cluster_list[[paste0("cluster", i)]] <- cluster_mean
}

expressionProfile_data <- as.data.frame(cluster_list)

# Create the design matrix from the `time_point` variable
des_mat <- model.matrix(~ 0 + MetaData$condition)
colnames(des_mat) <- levels(MetaData$condition)
# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(expressionProfile_data), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit2 <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit2, number = ncol(expressionProfile_data)) %>%
  tibble::rownames_to_column("cluster")

head(stats_df)

conditions <- c("S", "D", "SS", "DD", "SD", "DS")

# Generate all pairwise combinations
pairs <- combn(conditions, 2, simplify = FALSE)


contrast_defs <- sapply(pairs, function(pair) {
  paste0(pair[1], " - ", pair[2])
})

# Name each contrast like "S_vs_D"
names(contrast_defs) <- sapply(pairs, function(pair) {
  paste0(pair[1], "_vs_", pair[2])
})

contr.matrix <- makeContrasts(contrasts =contrast_defs,
  levels = colnames(des_mat))

fit2 <- contrasts.fit(fit, contr.matrix)
fit2 <- eBayes(fit2, robust = TRUE)

summary(decideTests(fit2, p.value = 0.1))
topTable(fit2, coef = "S_vs_D")

module_26_df <- expressionProfile_data %>%
  tibble::rownames_to_column("id") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(MetaData %>%
                      dplyr::select(id, condition),
                    by = c("id" = "id"))

ggplot(module_26_df,
       aes(
         x = condition,
         y = cluster1,
         color = condition)) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()



#### iterating through all clusters ####

# put the right number of clusters here
clusters <- paste0("cluster", 1:9)
conditions <- c("S", "D", "SS", "DD", "SD", "DS")
contrasts <- combn(conditions, 2, simplify = FALSE)
# Each element is like c("S", "D"), c("S", "SS"), etc.

# Initialize a list to hold per-cluster p-values
cluster_pval_lists <- vector("list", length(clusters))
names(cluster_pval_lists) <- clusters

# Loop over each cluster and get adjusted p-values for all contrasts
for (i in seq_along(clusters)) {
  clust <- clusters[i]
  pvals <- c()
  
  for (con in contrasts) {
    contrast_name <- paste0(con[1], " - ", con[2])  
    if (!(contrast_name %in% colnames(fit2$coefficients))) {
      contrast_name <- paste0(con[2], "_vs_", con[1])  # check reverse
    }
    
    tt <- topTable(fit2, coef = contrast_name, number = Inf, adjust.method = "BH")
    pval <- tt[clust, "adj.P.Val"]
    
    label <- paste(con, collapse = "-")
    pvals[label] <- pval
  }
  
  cluster_pval_lists[[clust]] <- pvals
}


# letters
cluster_letters <- lapply(cluster_pval_lists, function(pv) {
  multcompLetters(pv, threshold = 0.05)$Letters
})


eigengenes_long <- t(expressionProfile_data) %>%
  as.data.frame() %>%
  rownames_to_column("Cluster") %>%
  pivot_longer(-Cluster, names_to = "Sample", values_to = "Eigengene") %>%
  left_join(MetaData %>% select(Sample = id, Condition = condition), by = "Sample")

# individual plot
df1 <- eigengenes_long %>% filter(Cluster == "cluster1")

ggplot(df1, aes(x = Condition, y = Eigengene, color = Condition)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3)+
  #geom_jitter(width = 0.2) +
  theme_minimal() +
  geom_text(data = data.frame(
    Condition = names(cluster_letters$cluster1),
    label = cluster_letters$cluster1,
    y = tapply(df1$Eigengene, df1$Condition, max) + 0.05  # place letters above the boxes
  ),
  aes(x = Condition, y = y, label = label),
  vjust = 0
  ) +
  ggtitle("Cluster 1")


#### joint plots for all clusters ####
clusters_to_plot <- c("cluster10", "cluster9", "cluster8",  "cluster7", "cluster6", "cluster5", "cluster4", "cluster3","cluster2", "cluster1")
eigengenes_subset <- eigengenes_long %>% 
  filter(Cluster %in% clusters_to_plot)
# Convert cluster_letters list into long dataframe

annotation_letters <- bind_rows(
  lapply(names(cluster_letters), function(clust) {
    data.frame(
      Cluster = clust,
      Condition = names(cluster_letters[[clust]]),
      Letter = cluster_letters[[clust]]
    )
  }),
  .id = "id"
) %>% 
  filter(Cluster %in% clusters_to_plot)
# Calculate y position for placing letters above the boxplots
y_pos <- eigengenes_subset %>%
  group_by(Cluster, Condition) %>%
  summarise(y = max(Eigengene) + 0.05, .groups = "drop")

# Join with annotation letters
annotation_letters <- annotation_letters %>%
  left_join(y_pos, by = c("Cluster", "Condition"))
annotation_letters$Condition <- factor(annotation_letters$Condition, 
                                       levels = c("SD", "S", "SS", "D", "DD", "DS"))
eigengenes_subset$Condition <- factor(eigengenes_subset$Condition, 
                                      levels = c("SD", "S", "SS", "D", "DD", "DS"))


# annotation_letters$Condition <- factor(annotation_letters$Condition, 
#                                        levels = c("S", "SS", "SD", "D", "DD", "DS"))
# eigengenes_subset$Condition <- factor(eigengenes_subset$Condition, 
#                                        levels = c("S", "SS", "SD", "D", "DD", "DS"))


fill_colors <- c(
  "S"       = "#f75f55",
  "SS"    = "#f75f55",
  "SD"    = "#f75f55",     # color like 40
  "D"       = "#00A9FF",
  "DD"    = "#00A9FF",
  "DS"    = "#00A9FF"     # color like 10
)

fill_patterns <- c(
  "S"       = "none",
  "SS"    = "circle",
  "SD"    = "stripe",   # pattern like 10→10
  "D"       = "none",
  "DD"    = "stripe",
  "DS"    = "circle" # pattern like 40→40
)

bars <- ggplot(eigengenes_subset, aes(x = Condition, y = Eigengene)) +
  geom_boxplot_pattern(
    aes(fill = Condition, pattern = Condition),
    width = 0.1,
    #outlier.shape = NA,
    pattern_fill = "gray30",
    pattern_density = 0.3,
    pattern_spacing = 0.05,
    alpha = 1,
    lwd = 0.4
  ) +
  #geom_sina(aes(color = Condition), maxwidth = 0.3, alpha = 0.7) +
  geom_text(
    data = annotation_letters,
    aes(x = Condition, y = y, label = Letter),
    nudge_y = 0.05 * max(eigengenes_subset$Eigengene),
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Cluster, ncol = 3, scales = "free_y") +
  labs(
    title = "Eigengene expression across clusters",
    y = "Eigengene expression",
    x = "Condition"
  ) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = fill_patterns) +
  guides(fill = FALSE, pattern = FALSE, color = FALSE) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"), title = element_text(size = 11))
bars

#### drawing complex heatmap and bar plot together ####

# Heatmap object
ht_drawn <- draw(ht)

# ggplot barplot

# Convert ggplot to a grob
p_grob <- ggplotGrob(bars)

# Use grid layout
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(1, 1.2), "null"))))
# Draw heatmap in left column
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht_drawn, newpage = FALSE)
decorate_heatmap_body("Eigengene", {
  grid.text("Module\nclusters", x = unit(-12, "lines"), y = 1.2, rot = 0, gp = gpar(fontsize = 11, fontface = "italic"))
})
upViewport()

# Draw barplot in right column
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(p_grob)
upViewport()


###  Gene relationship to trait and important modules: Gene Significance and Module Membership

#Colors of the modules

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

SD <- as.data.frame(MetaData_bin)
names(SD)

geneTraitSignificance = as.data.frame(cor(datExpr, SD, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(SD), sep="")
names(GSPvalue) = paste("p.GS.", names(SD), sep="")

## Summary output of network analysis results
# Make a dataframe that connects traits, genes, and gene annotation

#Import annotation file.
annot_final <- read.csv("../Past_annot.csv", header = TRUE, sep = ",")
GO.annot <- subset(annot_final, select= c(SeqName, Length, GO_IDs)) #Select only relevant information

names(GO.annot)[1] <- "gene_id" #rename column
names(GO.annot)[2] <- "length"

#Match up genes in datExpr to annotation file
names(GO.annot)
probes = names(datExpr)
probes2annot = match(probes, GO.annot$gene_id)
# The following is the number of probes without annotation... Should return 0.
sum(is.na(probes2annot))

#Create the starting data frame
geneInfo0 = data.frame(gene_id = probes,
                       Annotation.GO.ID = GO.annot$GO_IDs[probes2annot],
                       moduleColor = paste0("ME",mergedColors),
                       geneTraitSignificance,
                       GSPvalue)

cond <- as.data.frame(as.numeric(MetaData$condition))
rownames(cond) <-  MetaData$id
#Order modules by their significance for time_point
modOrder = order(-abs(cor(MEs, SD$SD, use = "p")))

#Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

#Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.SD));
geneInfo = geneInfo0[geneOrder, ]
head(geneInfo)

#Add module cluster in geneInfo
geneInfo <- left_join(geneInfo, moduleCluster, by = "moduleColor")
dim(geneInfo)
head(geneInfo)

#Save geneInfo as a CSV
geneInfo$Annotation.GO.ID <- gsub(";NA", "", geneInfo$Annotation.GO.ID) #Remove NAs
geneInfo$Annotation.GO.ID <- gsub("NA", "", geneInfo$Annotation.GO.ID) #Remove NAs
write.csv(geneInfo, file = "rin/geneInfo.vst.outliers.csv")

#save.image(file="rin/070725.vst.RData") 
#load("rin/070725.vst.RData")


#### hub genes ####

# This gives correlation between each gene and each module eigengene
kME <- as.data.frame(cor(datExpr, MEs, use = "p"))
kME_pval <- as.data.frame(corPvalueStudent(as.matrix(kME), nSamples = nrow(datExpr)))

topHubGenes <- list()

allColors <- unique(mergedColors)
allColors <- allColors[allColors != "grey"]  # grey = unassigned

for (color in allColors) {
  # Gene indices for the current module
  inModule <- mergedColors == color
  
  # Column of the module eigengene in kME table
  column <- paste0("ME", color)
  
  # Subset kME and p-values for genes in the module
  kME_col <- kME[inModule, column]
  pval_col <- kME_pval[inModule, column]
  
  genes <- colnames(datExpr)[inModule]
  
  # Filter for significant genes only (e.g. p < 0.05)
  sig_idx <- which(pval_col < 0.05 & !is.na(pval_col) & !is.na(kME_col))
  
  if(length(sig_idx) == 0) {
    warning(paste0("No significant genes found in module ", color))
    next
  }
  
  kME_sig <- kME_col[sig_idx]
  genes_sig <- genes[sig_idx]
  
  # Order significant genes by absolute kME (can also use signed kME if preferred)
  topGenes <- genes_sig[order(abs(kME_sig), decreasing = TRUE)]
  
  # Save top 10 significant genes, or fewer if less than 10 significant genes
  topHubGenes[[color]] <- topGenes[1:min(200, length(topGenes))]
}

topHubGenes


hub_df <- do.call(rbind, lapply(names(topHubGenes), function(mod) {
  data.frame(module = mod,
             gene = topHubGenes[[mod]],
             rank = 1:length(topHubGenes[[mod]]))
}))
write.csv(hub_df, "rin/top_hub_genes_long_format.csv", row.names = FALSE)


hubs = chooseTopHubInEachModule(datExpr, mergedColors)


#### analysis of specific modules ####

# single module

# pairs and trios

# Long format
cor_long <- as.data.frame(moduleTraitCor) %>%
  rownames_to_column("module") %>%
  pivot_longer(-module, names_to = "trait", values_to = "correlation")

pval_long <- as.data.frame(moduleTraitPvalue) %>%
  rownames_to_column("module") %>%
  pivot_longer(-module, names_to = "trait", values_to = "pvalue")

# Join
trait_assoc <- left_join(cor_long, pval_long, by = c("module", "trait"))


# modules and traits correlation
trait_assoc %>%
  group_by(module) %>%
  summarise(
    n_strong = sum(abs(correlation) > 0.2 & pvalue < 0.05),
    traits = paste(trait[abs(correlation) > 0.2 & pvalue < 0.05], collapse = ", ")
  ) %>%
  filter(n_strong > 0) -> traitcor


write.csv(traitcor, "rin/module.trait.corr.csv", row.names = FALSE)

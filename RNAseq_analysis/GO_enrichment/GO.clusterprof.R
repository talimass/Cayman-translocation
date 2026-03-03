library(goseq)
library(tidyverse)
library(GSEABase)               #BiocManager::install("GSEABase")
library(data.table)
library(ggplot2)
library(cowplot)                #install.packages("cowplot")
library(patchwork)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(ontologyIndex)
library(GOSemSim)
library(simplifyEnrichment)
library(org.Hs.eg.db)
library(rlang)
library(readr)
library(stringr)
library(purrr)

# setting working directory 
setwd("/home/gospozha/haifa/cayman/rna/mapping/github")

# prepare term to gene table
eggNOG <- read.csv("../Past_annot.csv", header = TRUE, sep = ",") %>%
  dplyr::select(GO_IDs, SeqName) %>%
  dplyr::filter(GO_IDs!="NA;NA") %>%
  dplyr::filter(GO_IDs!=";;") %>%
  dplyr::filter(GO_IDs!="") %>%
  separate_rows(GO_IDs, sep = ";") %>%
  dplyr::select("term" = GO_IDs, "gene" = SeqName) %>%
  mutate(term = str_remove_all(term, '"')) %>%  # remove quotes
  filter(term != "") %>%
  distinct() %>%
  drop_na()

# prepare the term to name table
ontology <- get_ontology(file = "go.obo",
                         propagate_relationships = "is_a",
                         extract_tags = "everything",
                         merge_equivalent_terms = TRUE)
eggNOG_term <- eggNOG %>%
  mutate(name = ontology$name[term]) %>%
  select(c(term, name)) %>%
  distinct() %>%
  drop_na() %>%
  filter(!grepl("obsolete", name))

eggNOG <- eggNOG %>%
  filter(term %in% eggNOG_term$term)

# save the results
setwd("/home/gospozha/haifa/cayman/rna/mapping/github/rin/GO.new")
write_tsv(x = eggNOG, file = "term2gene_GO.tsv")
write_tsv(x = eggNOG_term, file = "term2name_GO.tsv")


countData  <- read.csv2('../../CountMatrix.csv', header=TRUE, row.names=1, sep=',', check.names = FALSE)
#gene count matrix
countData$geneID <- rownames(countData)

# filtering genes 
smallestGroupSize <- 4
keep <- rowSums(countData >= 10) >= smallestGroupSize
countData <- countData[keep,]

background_genes <- countData %>%
  dplyr::select("geneID") %>%
  unlist() %>%
  as.vector()

# perform ORA
term2gene <- read_tsv("term2gene_GO.tsv")
term2name <- read_tsv("term2name_GO.tsv")

#### DE genes ####
# List of comparisons
comparisons <- c("S.D", "SS.DD", "S.SS", "D.DD", "SS.SD", "SS.DS", "DD.DS", "DD.SD")

# Function to process one comparison
run_enrichment <- function(comp_name) {
  message("Processing: ", comp_name)
  
  # Load DEG list
  file_path <- paste0("../../GO/", comp_name, ".DE.genes.lfc1.csv")
  interesting_set <- read_csv(file_path, show_col_types = FALSE) %>%
    dplyr::select("geneID" = 1) %>%
    pull()
  
  # Run enrichment
  enrichment <- enricher(interesting_set,
                         TERM2GENE = term2gene,
                         TERM2NAME = term2name,
                         pvalueCutoff = 0.05,
                         universe = background_genes,
                         pAdjustMethod = "fdr",
                         qvalueCutoff = 0.2)
  
  # Save enrichment results
  write_csv(enrichment@result, paste0("enrichment_results.", comp_name, ".lfc1.csv"))

# Run for all comparisons
lapply(comparisons, run_enrichment)


#### bar plot ####
deg_results <- list(
  SvD   = read_csv("../GO/S.D.DE.genes.csv"),
  SSvDD = read_csv("../GO/SS.DD.DE.genes.csv"),
  SvSS  = read_csv("../GO/S.SS.DE.genes.csv"),
  DvDD  = read_csv("../GO/D.DD.DE.genes.csv"),
  SSvSD = read_csv("../GO/SS.SD.DE.genes.csv"),
  SSvDS = read_csv("../GO/SS.DS.DE.genes.csv"),
  DDvDS = read_csv("../GO/DD.DS.DE.genes.csv"),
  DDvSD = read_csv("../GO/DD.SD.DE.genes.csv")
)

# Enrichment results
enrich_results <- list(
  SvD   = read_csv("enrichment_results.S.D.csv"),
  SSvDD = read_csv("enrichment_results.SS.DD.csv"),
  SvSS  = read_csv("enrichment_results.S.SS.csv"),
  DvDD  = read_csv("enrichment_results.D.DD.csv"),
  SSvSD = read_csv("enrichment_results.SS.SD.csv"),
  SSvDS = read_csv("enrichment_results.SS.DS.csv"),
  DDvDS = read_csv("enrichment_results.DD.DS.csv"),
  DDvSD = read_csv("enrichment_results.DD.SD.csv")
)


# Function to compute mean logFC per GO term
compute_mean_logfc <- function(enrich_df, deg_df) {
  enrich_df %>%
    select(ID, Description, geneID, p.adjust) %>%
    separate_rows(geneID, sep = "/") %>%
    left_join(deg_df, by = c("geneID" = "...1")) %>%
    group_by(ID, Description, p.adjust) %>%
    dplyr::summarize(mean_logFC = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")
}

# Apply to all comparisons
mean_logfc_list <- map2(enrich_results, deg_results, compute_mean_logfc)

# Add comparison name
mean_logfc_named <- map2_dfr(mean_logfc_list, names(mean_logfc_list), ~mutate(.x, contrast = .y))

# Reverse DDvDS to DSvDD
mean_logfc_named <- mean_logfc_named %>%
  mutate(
    contrast = ifelse(contrast == "DDvDS", "DSvDD", contrast),
    mean_logFC = ifelse(contrast == "DSvDD", -mean_logFC, mean_logFC)
  )%>%
  mutate(
    contrast = ifelse(contrast == "DDvSD", "SDvDD", contrast),
    mean_logFC = ifelse(contrast == "SDvDD", -mean_logFC, mean_logFC)
  )

# You can optionally filter top N GO terms across all contrasts
top_terms <- mean_logfc_named %>%
  mutate(z_logFC = scale(mean_logFC)[,1]) %>%
  dplyr::filter(p.adjust<0.05)%>%
  mutate(Description_ordered = factor((Description), levels = rev
                                      (unique(Description))))

# the same colors
# Define your contrasts in desired order
contrast_levels <- (c("DSvDD","SDvDD","SSvDS","SSvSD","SSvDD","SvD","DvDD","SvSS"))

top_terms$contrast <- factor(top_terms$contrast, levels = contrast_levels)

# Define colors and shapes
contrast_colors <- c(
  "SvSS" = "#1f78b4",
  "DvDD" = "#33a02c",
  "SvD" = "#e31a1c",
  "SSvDD" = "#ff7f00",
  "SSvSD" = "#6a3d9a",
  "DSvDD" = "#b15928",
  "SSvDS" = "#a6cee3",
  "SDvDD" = "#fb9a99"
)

contrast_shapes <- c(
  "SvSS" = 21,
  "DvDD" = 22,
  "SvD" = 23,
  "SSvDD" = 24,
  "SSvSD" = 25,
  "DSvDD" = 8,
  "SSvDS" = 3,
  "SDvDD" = 4
)


# Plot barplot
# either z score of logfc or pure logfc
clust.barplot <- ggplot(
  top_terms,
  aes(x = Description_ordered,
    y = mean_logFC,
    fill = contrast,
    shape = contrast)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.6) +
  geom_point(size = 3, position = position_dodge(width = 0.8), color = "black") +
  coord_flip() +
  #facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
  scale_y_continuous(
    name = "Mean log2FC\n(Down in S  ←  0  →  Up in S)",
    breaks = scales::breaks_width(5)
  ) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 35), name = "GO term")+
  scale_fill_manual(values = contrast_colors) +
  scale_shape_manual(values = contrast_shapes) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) +
  #labs(title = "GO term enrichment") +
  theme(
    plot.title = element_text(size = 13, margin = margin(t = 9, b = 6)),
    plot.margin = margin(10, 10, 8, 8),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    legend.position = "right"
  )



clust.barplot

ggsave("../../fin/clust.barplot.not.ordered.jpg", clust.barplot, width = 12, height = 7.5)


#### separate up and down ####
run_enrichment_split <- function(comp_name) {
  message("Processing: ", comp_name)
  
  file_path <- paste0("../../GO/", comp_name, ".DE.genes.lfc1.csv")
  deg <- read_csv(file_path, show_col_types = FALSE)
  
  # --- adjust these column names if yours differ ---
  # gene column
  if (!"geneID" %in% names(deg)) names(deg)[1] <- "geneID"
  # log2FC column (common names)
  if (!"log2FoldChange" %in% names(deg)) {
    stop("No log2FoldChange column found in: ", file_path)
  }
  
  up_genes <- deg %>% filter(log2FoldChange > 0) %>% pull(geneID) %>% unique()
  down_genes <- deg %>% filter(log2FoldChange < 0) %>% pull(geneID) %>% unique()
  
  enrich_one <- function(gene_vec) {
    if (length(gene_vec) < 5) return(NULL)  # too few genes to be meaningful
    enricher(gene_vec,
             TERM2GENE = term2gene,
             TERM2NAME = term2name,
             pvalueCutoff = 0.05,
             universe = background_genes,
             pAdjustMethod = "fdr",
             qvalueCutoff = 0.2)
  }
  
  enrich_up <- enrich_one(up_genes)
  enrich_down <- enrich_one(down_genes)
  
  # write results (even if empty, write a blank tibble for consistency)
  up_out <- if (!is.null(enrich_up)) enrich_up@result else tibble()
  down_out <- if (!is.null(enrich_down)) enrich_down@result else tibble()
  
  write_csv(up_out,   paste0("enrichment_results.", comp_name, ".UP.lfc1.csv"))
  write_csv(down_out, paste0("enrichment_results.", comp_name, ".DOWN.lfc1.csv"))
  
  invisible(list(UP = enrich_up, DOWN = enrich_down))
}

lapply(comparisons, run_enrichment_split)

#### bar plot for separated up and down ####
deg_results <- list(
  SvD   = read_csv("../GO/S.D.DE.genes.csv"),
  SSvDD = read_csv("../GO/SS.DD.DE.genes.csv"),
  SvSS  = read_csv("../GO/S.SS.DE.genes.csv"),
  DvDD  = read_csv("../GO/D.DD.DE.genes.csv"),
  SSvSD = read_csv("../GO/SS.SD.DE.genes.csv"),
  SSvDS = read_csv("../GO/SS.DS.DE.genes.csv"),
  DDvDS = read_csv("../GO/DD.DS.DE.genes.csv"),
  DDvSD = read_csv("../GO/DD.SD.DE.genes.csv")
)


enrich_files <- list.files(
  "./UPandDOWN/altered",
  pattern = "^enrichment_results\\..*\\.(UP|DOWN)\\.lfc1\\.csv$",
  full.names = TRUE
)

enrich_results <- map(enrich_files, read_csv, show_col_types = FALSE)

names(enrich_results) <- enrich_files %>%
  basename() %>%
  str_remove("^enrichment_results\\.") %>%   # remove prefix
  str_remove("\\.lfc1\\.csv$")               # remove suffix

compute_mean_logfc <- function(enrich_df, deg_df) {
  enrich_df %>%
    select(ID, Description, geneID, p.adjust) %>%
    separate_rows(geneID, sep = "/") %>%
    left_join(deg_df, by = c("geneID" = "...1")) %>%
    group_by(ID, Description, p.adjust) %>%
    summarise(
      mean_logFC = mean(log2FoldChange, na.rm = TRUE),
      .groups = "drop"
    )
}

mean_logfc_named <- imap_dfr(enrich_results, function(enrich_df, nm) {
  
  # nm example: "S.D.UP.lfc1" or "S.D.UP" depending on how you named it
  m <- str_match(nm, "^([^.]+\\.[^.]+)\\.(UP|DOWN)$")
  if (is.na(m[1,1])) {
    warning("Skipping (unexpected name): ", nm)
    return(tibble())
  }
  
  contrast_raw <- m[1,2]   # "S.D"
  direction    <- m[1,3]   # "UP" or "DOWN"
  
  # Convert "S.D" -> "SvD" etc to match your deg_results names
  contrast_key <- str_replace(contrast_raw, "\\.", "v")
  
  deg_df <- deg_results[[contrast_key]]
  if (is.null(deg_df)) {
    warning("No DEG table found for contrast_key = ", contrast_key, " (from file name: ", nm, ")")
    return(tibble())
  }
  
  out <- compute_mean_logfc(enrich_df, deg_df)
  
  # DOWN gene set corresponds to negative log2FC in DE table -> flip for plotting consistency
  #if (direction == "DOWN") out$mean_logFC <- -out$mean_logFC
  
  out %>% mutate(contrast = contrast_key, direction = direction)
})


mean_logfc_named <- mean_logfc_named %>%
  mutate(
    contrast = ifelse(contrast == "DDvDS", "DSvDD", contrast),
    mean_logFC = ifelse(contrast == "DSvDD", -mean_logFC, mean_logFC)
  ) %>%
  mutate(
    contrast = ifelse(contrast == "DDvSD", "SDvDD", contrast),
    mean_logFC = ifelse(contrast == "SDvDD", -mean_logFC, mean_logFC)
  )


top_terms <- mean_logfc_named %>%
  filter(p.adjust < 0.05) %>%
  mutate(
    z_logFC = scale(mean_logFC)[,1],
    Description_ordered = factor(
      Description,
      levels = rev(unique(Description))
    )
  )


# the same colors
# Define your contrasts in desired order
contrast_levels <- (c("DSvDD","SDvDD","SSvDS","SSvSD","SSvDD","SvD","DvDD","SvSS"))

top_terms$contrast <- factor(top_terms$contrast, levels = contrast_levels)

# Define colors and shapes
contrast_colors <- c(
  "SvSS" = "#1f78b4",
  "DvDD" = "#33a02c",
  "SvD" = "#e31a1c",
  "SSvDD" = "#ff7f00",
  "SSvSD" = "#6a3d9a",
  "DSvDD" = "#b15928",
  "SSvDS" = "#a6cee3",
  "SDvDD" = "#fb9a99"
)

contrast_shapes <- c(
  "SvSS" = 21,
  "DvDD" = 22,
  "SvD" = 23,
  "SSvDD" = 24,
  "SSvSD" = 25,
  "DSvDD" = 8,
  "SSvDS" = 3,
  "SDvDD" = 4
)


# Plot barplot
# either z score of logfc or pure logfc

clust.barplot <- ggplot(
  top_terms,
  aes(x = Description_ordered,
      y = mean_logFC,
      fill = contrast,
      shape = contrast)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, alpha = 0.6) +
  geom_point(size = 2, position = position_dodge(width = 0.8), color = "black") +
  coord_flip() +
  #facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
  scale_y_continuous(
    name = "Mean log2FC\n(Down in S  ←  0  →  Up in S)",
    breaks = scales::breaks_width(5)
  ) +
  scale_x_discrete(
    labels = function(x) stringr::str_trunc(x, width = 50),
    name = "GO term"
  )+
  scale_fill_manual(values = contrast_colors) +
  scale_shape_manual(values = contrast_shapes) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) +
  labs(title = "GO term enrichment") +
  theme(
    plot.title = element_text(size = 13, margin = margin(t = 9, b = 6)),
    plot.margin = margin(10, 10, 8, 8),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 8),
    legend.position = "right"
  )



clust.barplot
ggsave("../../fin/clust.barplot.separated.pdf", clust.barplot, width = 5, height = 12)


#### making figure 5 ####

deg.bar.diverging <- readRDS("~/haifa/cayman/rna/mapping/github/fin/deg.bar.RDS")


combined_go <- (deg.bar.diverging | clust.barplot) +
  plot_layout(widths = c(1.2, 0.8), heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A', tag_prefix = '', tag_suffix = '')

#theme(plot.margin = margin(10, 10, 10, 10))
combined_go
# display it
ggsave("~/haifa/cayman/rna/mapping/github/fin/combined.fig5.jpg", combined_go, width = 12, height =7)
ggsave("~/haifa/cayman/rna/mapping/github/fin/combined.fig5.pdf", combined_go, width = 12, height =7)


#### combined fig 5 and 6 ####
biomin.bar <- readRDS("~/haifa/cayman/rna/mapping/github/fin/biomin.bar.RDS")
biomin.bar <- biomin.bar +
  labs(tag = "C.") +
  theme(
    plot.tag = element_text(size = 13, face = "plain"),
    plot.tag.position = c(0, 0.969)
  )



combined_go <-
  (deg.bar.diverging | ( clust.barplot /biomin.bar )) +
  plot_layout(
    widths  = c(0.4, 0.6),  # controls C width
    heights = c(1, 1))+
  plot_annotation(tag_levels = 'A', tag_prefix = '', tag_suffix = '.')

combined_go

ggsave("~/haifa/cayman/rna/mapping/github/fin/combined.fig56.jpg", combined_go, width = 12, height =10)
ggsave("~/haifa/cayman/rna/mapping/github/fin/combined.fig56.pdf", combined_go, width = 12, height =7)


#### WGCNA Modules ####
setwd("/home/gospozha/haifa/cayman/rna/mapping/github/rin/GO.new")
# Read your data
# Get module name from moduleColor, e.g. "MEblack" -> "black"
# Create the correct p-value column name for each gene
# Now use rowwise() to pull out the correct pval for each gene
cluster.wgcna <- read.csv( "../geneInfo.vst.outliers.csv", header = TRUE, row.names = 1) %>%
  mutate(module = gsub("^ME", "", moduleColor)) %>%
  mutate(pval_col = paste0("p.MM.", module), mm_col = (paste0("MM.", module))) %>%
  rowwise() %>%
  mutate(module_pval = get(pval_col), module_mm = get(mm_col)) %>%
  ungroup()

# Filter for significant genes (module p < 0.05)
signif_genes <- cluster.wgcna %>%
  filter(module_pval < 0.05, abs(module_mm) > 0.8)

# add module GS > 0.2, p val < 0.05 at least in one column
# 19,543 rlog
# 8169
# 6636 vst outliers

# filtering by GS
# Define all trait names used in your dataset
traits <- c("S", "D", "SD", "DD", "SS", "DS")

# Build a filter expression:
gs_filter <- purrr::map2_chr(
  traits,
  paste0("p.GS.", traits),
  ~ paste0("(abs(GS.", .x, ") > 0.2 & ", .y, " < 0.05)")
) %>%
  paste(collapse = " | ")  # join with OR

# Combine the module and GS filters using dplyr + rlang

signif_genes <- cluster.wgcna %>%
  filter(module_pval < 0.05, abs(module_mm) > 0.8) %>%
  filter(!!parse_expr(gs_filter))

#6,394 genes

clustered_signif_genes <- signif_genes %>%
  group_by(moduleCluster) %>%
  summarise(genes = list(gene_id))

# List of comparisons
comparisons <- c(1:9)

# Function to process one comparison
run_enrichment <- function(comp_name) {
  message("Processing: ", comp_name)

  interesting_set <- clustered_signif_genes %>%
    dplyr::filter(moduleCluster == comp_name) %>%
    dplyr::select(genes)%>%
    unlist()%>%
    as.factor()
  
  # Run enrichment
  enrichment <- enricher(interesting_set,
                         TERM2GENE = term2gene,
                         TERM2NAME = term2name,
                         pvalueCutoff = 0.05,
                         universe = background_genes,
                         pAdjustMethod = "fdr",
                         qvalueCutoff = 0.2)
  enrichment@result <- enrichment@result %>% 
    dplyr::filter(Count >= 2)
  write_csv(enrichment@result, paste0("enrichment_results.cluster", comp_name, ".vst.GS2.csv"))  
  # Plot dotplot if any significant terms
  
  if (any(enrichment@result$p.adjust <= 0.05)) {
    # Save enrichment results

    p <- clusterProfiler::dotplot(enrichment,
                 x = "geneRatio",
                 color = "p.adjust",
                 orderBy = "x",
                 showCategory = 100,
                 font.size = 8) +
      ggtitle(paste0("Cluster ", comp_name))
    
    ggsave(paste0("enrichment_dotplot.cluster", comp_name, ".vst.outliers2.jpg"), p, width = 4, height = 6)
  }
}

# Run for all comparisons
lapply(comparisons, run_enrichment)

# more beautiful plot for each module

# List of comparisons
comparisons <- c(4)

# Function to process one comparison
run_enrichment <- function(comp_name) {
  message("Processing: ", comp_name)
  
  interesting_set <- clustered_signif_genes %>%
    dplyr::filter(moduleCluster == comp_name) %>%
    dplyr::select(genes)%>%
    unlist()%>%
    as.factor()
  
  # Run enrichment
  enrichment <- enricher(interesting_set,
                         TERM2GENE = term2gene,
                         TERM2NAME = term2name,
                         pvalueCutoff = 0.05,
                         universe = background_genes,
                         pAdjustMethod = "fdr",
                         qvalueCutoff = 0.2)
  enrichment@result <- enrichment@result %>% 
    dplyr::filter(Count >= 2)
  write_csv(enrichment@result, paste0("enrichment_results.cluster", comp_name, ".vst.GS2.csv"))  
  # Plot dotplot if any significant terms
  
  if (any(enrichment@result$p.adjust <= 0.05)) {
    # Save enrichment results
    
    p <- clusterProfiler::dotplot(enrichment,
                                  x = "geneRatio",
                                  color = "p.adjust",
                                  orderBy = "x",
                                  showCategory = 100,
                                  font.size = 7.5) +
      ggtitle(paste0("WGCNA: high correlation with deep"))+
      theme(
        axis.text = element_text(size = 8),          # axis tick labels
        axis.title = element_text(size = 12),         # axis titles
        legend.title = element_text(size = 11),       # legend title
        legend.text = element_text(size = 9)          # legend labels
      )
    
    ggsave(paste0("enrichment_dotplot.cluster", comp_name, ".new.jpg"), p, width = 4, height = 9)
  }
}

# Run for all comparisons
lapply(comparisons, run_enrichment)



#### heatmap wgcna for modules ####

cluster_gene_lists <- signif_genes %>%
  #filter(!moduleCluster %in% c(8)) %>%  # remove clusters you don't want to show
  group_by(moduleCluster) %>%
  summarise(genes = list(gene_id), .groups = "drop") %>%
  deframe()


# Enrichment results
enrich_results <- list(
  cluster1 = read_csv("enrichment_results.cluster1.vst.GS.csv"),
  cluster2 = read_csv("enrichment_results.cluster2.vst.GS.csv"),
  cluster3 = read_csv("enrichment_results.cluster3.vst.GS.csv"),
  cluster4 = read_csv("enrichment_results.cluster4.vst.GS.csv"),
  cluster5 = read_csv("enrichment_results.cluster5.vst.GS.csv"),
  cluster6 = read_csv("enrichment_results.cluster6.vst.GS.csv"),
  cluster7 = read_csv("enrichment_results.cluster7.vst.GS.csv"),
  cluster8 = read_csv("enrichment_results.cluster8.vst.GS.csv"),
  cluster9 = read_csv("enrichment_results.cluster9.vst.GS.csv")
)

# Enrichment results (with removed mammal-derived factors)
enrich_results <- list(
  cluster1 = read_csv("enrichment_results.cluster1.vst.GS.csv"),
  cluster2 = read_csv("enrichment_results.cluster2.vst.cor.GS.csv"),
  cluster3 = read_csv("enrichment_results.cluster3.vst.GS.csv"),
  cluster4 = read_csv("enrichment_results.cluster4.vst.GS.csv"),
  cluster5 = read_csv("enrichment_results.cluster5.vst.cor.GS.csv"),
  cluster6 = read_csv("enrichment_results.cluster6.vst.cor.GS.csv"),
  cluster7 = read_csv("enrichment_results.cluster7.vst.GS.csv"),
  #cluster8 = read_csv("enrichment_results.cluster8.vst.GS.csv"),
  cluster9 = read_csv("enrichment_results.cluster9.vst.cor.GS.csv")
)

# Enrichment results 26 sep
enrich_results <- list(
  cluster1 = read_csv("enrichment_results.cluster1.vst.GS2.csv"),
  cluster2 = read_csv("enrichment_results.cluster2.vst.cor.GS2.csv"),
  cluster3 = read_csv("enrichment_results.cluster3.vst.cor.GS2.csv"),
  cluster4 = read_csv("enrichment_results.cluster4.vst.cor.GS2.csv"),
  cluster5 = read_csv("enrichment_results.cluster5.vst.GS2.csv"),
  cluster6 = read_csv("enrichment_results.cluster6.vst.cor.GS2.csv"),
  cluster7 = read_csv("enrichment_results.cluster7.vst.cor.GS2.csv")
  #cluster8 = read_csv("enrichment_results.cluster8.vst.GS2.csv"),
  #cluster9 = read_csv("enrichment_results.cluster9.vst.GS2.csv")
)


all_enrich <- imap_dfr(enrich_results, function(df, cluster) {
  df %>%
    filter(p.adjust < 0.05) %>%
    mutate(cluster = cluster) %>%
    select(ID, Description, cluster, Count, p.adjust)
})

all_enrich <- all_enrich %>%
  mutate(logp = -log10(p.adjust)) %>%
  mutate(Description = str_trunc(Description, 60))  # limit term text

  
plot_data <- all_enrich %>%
  mutate(Description = factor(Description, levels = unique(Description)))

heatmap_clusters2 <- ggplot(plot_data, aes(x = cluster, y = Description, fill = logp)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(x = "Cluster", y = "GO Term", fill = "-log10(p.adj)", title = "Enriched GO Terms per WGCNA cluster") +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 50, hjust = 1, size = 9), plot.title = element_text(hjust=0.8, size = 13))

heatmap_clusters2
ggsave("../../fin/heatmap_clusters26sep.cor.jpg", heatmap_clusters2, width = 6, height = 12)


save.image(file="260925.GO.RData") 
#load("rin/070725.vst.RData")

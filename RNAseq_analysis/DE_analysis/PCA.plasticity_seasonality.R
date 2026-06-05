library(vegan)
library(ggplot2)
library(dplyr)
library(rcompanion)
library(performance)
library(emmeans)
library(DESeq2)
library(limma)
library(vegan)
library(pairwiseAdonis)
library(rcompanion)
library(multcompView)

# setting working directory 
setwd("/home/gospozha/haifa/cayman/rna/mapping/github/plasticity/")

#### PCA on translocation ####

# reading count matrix from a file
countData  <- read.csv2('../CountMatrix.csv', header=TRUE, row.names=1, sep=',', check.names = F)
# reading metadata file
MetaData <- read.csv2('../Metadata.treatment.csv', header=TRUE, sep=',')

MetaData$condition <- factor(MetaData$condition, levels = c("SS", "SD", "DD", "DS"))
MetaData$site <- factor(MetaData$site, levels = c("MF","CC"))
MetaData$RIN2 <- as.factor(MetaData$RIN2)
MetaData$rRNA <- as.factor(MetaData$rRNA)
MetaData$origin_depth <- as.factor(MetaData$origin_depth)
MetaData$transplant_depth <- as.factor(MetaData$transplant_depth)

# sample names in both objects
samples_meta  <- MetaData$id
samples_count <- colnames(countData)
# find common samples
common_samples <- intersect(samples_meta, samples_count)
# subset and reorder count matrix
countData<- countData[, common_samples]
# reorder metadata to match countData
MetaData <- MetaData[match(common_samples, MetaData$id), ]
# must be TRUE
all(colnames(countData) == MetaData$id)

# filtering for low counts
keep <- rowSums(countData >= 10) >= 4

countData_filt <- countData[keep, ]

cat("Genes retained:", nrow(countData_filt), "\n")

# vst transformation

dds <- DESeqDataSetFromMatrix(
  countData = countData_filt,
  colData = MetaData,
  design = ~ site  + rRNA + RIN2 + condition 
)

dds <- estimateSizeFactors(dds)

vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# remove batch effect
mm <- model.matrix(~condition + site  + rRNA + RIN2 , colData(vsd))
treatment.design <- mm[,1:4]
batch.design <- mm[,-(1:4)]
mat <- limma::removeBatchEffect(vsd_mat, covariates=batch.design, design=treatment.design)
vsd_filtered <- mat

# PCA
pca_input <- t(as.matrix(vsd_filtered))
pca <- prcomp(pca_input, scale. = FALSE) 

percent_var <- pca$sdev^2 / sum(pca$sdev^2)
percent_var_percent <- percent_var * 100

# Scree table
pca_variance_table <- tibble(
  PC = paste0("PC", seq_along(percent_var)),
  variance_explained = percent_var,
  percent_variance = percent_var_percent,
  cumulative_percent = cumsum(percent_var_percent)
)

head(pca_variance_table, 10)
#n_pcs <- which(cumsum(percent_var_percent) >= 70)[1]
n_pcs <- 2
pc_cols <- paste0("PC", 1:n_pcs)
pc_weights <- percent_var[1:n_pcs]

# Normalize weights so they sum to 1 among selected PCs
pc_weights <- pc_weights / sum(pc_weights)

pc_weights


pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("id") %>%
  left_join(MetaData, by = "id")

# Keep a clean group column
pca_df <- pca_df %>%
  mutate(
    group = condition
  )

head(pca_df)

#### PERMANOVA ####
# PERMANOVA on the weighted distances

# PERMANOVA on the complete batch-corrected transcriptome profile
# Transpose 
global_dist <- dist(t(vsd_filtered)) 

# PERMANOVA
permanova_res <- adonis2(global_dist ~ origin_depth * transplant_depth, 
                         data = MetaData, 
                         permutations = 999)
print(permanova_res)

# Marginal testing 
permanova_terms <- adonis2(global_dist ~ origin_depth * transplant_depth, 
                           data = MetaData, 
                           permutations = 999,
                           by = "terms") 
print(permanova_terms)


pairwise_res <- pairwise.adonis2(global_dist ~ condition, data = MetaData, permutations = 999)

# Print pairwise comparisons
print(pairwise_res)


# multivariate dispersion based on experimental groups
disp_res <- betadisper(global_dist, group = MetaData$condition)

# test
anova(disp_res)
permutest(disp_res, permutations = 999)
# betadisper not significant - PERMANOVA is trustworthy 

# permanova for plot text
p_orig <- permanova_terms["origin_depth", "Pr(>F)"]
p_trans <- permanova_terms["transplant_depth", "Pr(>F)"]
p_text <- paste0("P_Origin = ", p_orig,"\nP_Transplant = ", p_trans) 

#### PCA distance to centroids ####
centroids_2d <- pca_df %>%
  group_by(condition) %>%
  summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  )

arrow_df <- tibble(
  shift = c("SS to SD", "DD to DS"),
  x = c(
    centroids_2d$PC1[centroids_2d$condition == "SS"],
    centroids_2d$PC1[centroids_2d$condition == "DD"]
  ),
  y = c(
    centroids_2d$PC2[centroids_2d$condition == "SS"],
    centroids_2d$PC2[centroids_2d$condition == "DD"]
  ),
  xend = c(
    centroids_2d$PC1[centroids_2d$condition == "SD"],
    centroids_2d$PC1[centroids_2d$condition == "DS"]
  ),
  yend = c(
    centroids_2d$PC2[centroids_2d$condition == "SD"],
    centroids_2d$PC2[centroids_2d$condition == "DS"]
  )
)

# for plotting - inches to pt
pt_to_mm <- function(pt) pt / 2.845

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 1.5) +
  #stat_ellipse(aes(group = condition), linewidth = 0.5, alpha = 0.35, show.legend = FALSE) +
  geom_point(
    data = centroids_2d,
    aes(x = PC1, y = PC2, color = condition),
    inherit.aes = FALSE,
    shape = 18,
    size = 2.3
  ) +
  geom_segment(
    data = arrow_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.1, "cm")),
    linewidth = 0.4,
    color = "black"
  ) +
  geom_text_repel(
    data = centroids_2d,
    aes(x = PC1, y = PC2, label = condition, color = condition),
    inherit.aes = FALSE,
    size = pt_to_mm(7),
    show.legend = FALSE
  ) +
  theme_classic() +
  labs(
    x = paste0("PC1: ", round(percent_var_percent[1], 1), "% variance"),
    y = paste0("PC2: ", round(percent_var_percent[2], 1), "% variance")
    #title = "PCA of VST-normalized host gene expression"
  ) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )+
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = p_text,
    hjust = -0.05,
    vjust = 1.2,
    size = pt_to_mm(6)
  )
p_pca

####. Variance-weighted distance function ####
weighted_pc_distance <- function(sample_scores, centroid_scores, weights) {
  sqrt(sum(weights * (sample_scores - centroid_scores)^2))
}

# Calculate SS and DD native/control centroids

centroids <- pca_df %>%
  filter(condition %in% c("SS", "DD")) %>%
  group_by(condition) %>%
  summarise(
    across(all_of(pc_cols), mean, na.rm = TRUE),
    .groups = "drop"
  )

SS_centroid <- centroids %>%
  filter(condition == "SS") %>%
  dplyr::select(all_of(pc_cols)) %>%
  as.numeric()

DD_centroid <- centroids %>%
  filter(condition == "DD") %>%
  dplyr::select(all_of(pc_cols)) %>%
  as.numeric()

# Calculate distances and indices

distance_df <- pca_df %>%
  rowwise() %>%
  mutate(
    dist_to_SS = weighted_pc_distance(
      c_across(all_of(pc_cols)),
      SS_centroid,
      pc_weights
    ),
    dist_to_DD = weighted_pc_distance(
      c_across(all_of(pc_cols)),
      DD_centroid,
      pc_weights
    )
  ) %>%
  ungroup() %>%
  mutate(
    origin_centroid = case_when(
      origin_depth == "S" ~ "SS",
      origin_depth == "D" ~ "DD",
      TRUE ~ NA_character_
    ),
    destination_centroid = case_when(
      transplant_depth == "S" ~ "SS",
      transplant_depth == "D" ~ "DD",
      TRUE ~ NA_character_
    ),
    distance_from_origin = case_when(
      origin_centroid == "SS" ~ dist_to_SS,
      origin_centroid == "DD" ~ dist_to_DD,
      TRUE ~ NA_real_
    ),
    distance_to_destination = case_when(
      destination_centroid == "SS" ~ dist_to_SS,
      destination_centroid == "DD" ~ dist_to_DD,
      TRUE ~ NA_real_
    ),
    deep_like_index = dist_to_SS / (dist_to_SS + dist_to_DD)
  )


#### plasticity index ####
pca_plasticity_df <- distance_df %>%
  mutate(
    condition = factor(condition, levels = c("SS", "SD", "DD", "DS")),
    Plasticity = distance_from_origin
  )

# Check statistics
pca_plasticity_df %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean = mean(Plasticity, na.rm = TRUE),
    sd = sd(Plasticity, na.rm = TRUE),
    .groups = "drop"
  )


fit_aov <- aov(Plasticity ~ condition, data = pca_plasticity_df)

summary(fit_aov)
TukeyHSD(fit_aov)

# Assumption checks
shapiro.test(residuals(fit_aov))
car::leveneTest(Plasticity ~ condition, data = pca_plasticity_df)
# failed, repeated with kruskall

kruskal.test(Plasticity ~ condition, data = pca_plasticity_df)

pairwise.wilcox.test(
  pca_plasticity_df$Plasticity,
  pca_plasticity_df$condition,
  p.adjust.method = "BH"
)

# letters
pairwise_res <- pairwise.wilcox.test(
  pca_plasticity_df$Plasticity,
  pca_plasticity_df$condition,
  p.adjust.method = "BH"
)

# Convert pairwise p-value matrix to letters
pw_p <- pairwise_res$p.value

# multcompLetters needs named vector
pw_vec <- pw_p[!is.na(pw_p)]
names(pw_vec) <- apply(which(!is.na(pw_p), arr.ind = TRUE), 1, function(i) {
  paste(rownames(pw_p)[i[1]], colnames(pw_p)[i[2]], sep = "-")
})

cld <- multcompLetters(pw_vec)

letter_df <- data.frame(
  condition = names(cld$Letters),
  Letter = cld$Letters
)


plas_stats <- pca_plasticity_df %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean_p = mean(Plasticity, na.rm = TRUE),
    sd_p = sd(Plasticity, na.rm = TRUE),
    se_p = sd_p / sqrt(n),
    ci_val = 1.96 * se_p,
    lower = mean_p - ci_val,
    upper = mean_p + ci_val,
    max_group = max(Plasticity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(letter_df, by = "condition")

# Global y-position for significance letters
y_label <- max(pca_plasticity_df$Plasticity, na.rm = TRUE) * 1.12


p_plas <- ggplot(
  pca_plasticity_df,
  aes(x = condition, y = Plasticity, color = condition)
) +
  # Raw data points
  geom_jitter(
    width = 0.15,
    size = 1.5,
    #alpha = 0.65,
   # shape = 1,
   # stroke = 1
  ) +
  # Mean point
  geom_point(
    data = plas_stats,
    aes(x = condition, y = mean_p),
    color = "black",
    size = 1.5,
    inherit.aes = FALSE
  ) +
  
  # 95% CI
  geom_errorbar(
    data = plas_stats,
    aes(x = condition, ymin = lower, ymax = upper),
    color = "black",
    width = 0,
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  
  # Significance letters
  geom_text(
    data = plas_stats,
    aes(x = condition, y = y_label, label = Letter),
    color = "black",
    size = pt_to_mm(6),
    inherit.aes = FALSE
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  theme_classic() +
  labs(
    x = "Treatment",
    y = "Distance from origin"
   # title = "Weighted PCA transcriptomic plasticity"
  ) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    #axis.text.x = element_text(face = "bold"),
    legend.position = "none"
  )

p_plas


combined <- p_pca + p_plas
combined

#### deep index ####

p_deepness <- ggplot(distance_df, aes(x = condition, y = deep_like_index, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(
    x = "Treatment",
    y = "Deep-like transcriptomic index",
    title = "Transcriptomic similarity to shallow and deep controls"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

p_deepness

# Deep-like index across all four treatments
model_deepness <- lm(deep_like_index ~ condition, data = distance_df)
anova(model_deepness)
TukeyHSD(aov(model_deepness))

# Assumption checks
shapiro.test(residuals(model_deepness))


#### PCA on season ####

# reading count matrix from a file
countData  <- read.csv2('../CountMatrix.csv', header=TRUE, row.names=1, sep=',', check.names = F)
# reading metadata file
MetaData <- read.csv2('../Metadata.control.csv', header=TRUE, sep=',')

MetaData$condition <- factor(MetaData$condition, levels = c("S", "SS", "D", "DD"))
MetaData$site <- factor(MetaData$site, levels = c("MF","CC"))
MetaData$RIN2 <- as.factor(MetaData$RIN2)
MetaData$rRNA <- as.factor(MetaData$rRNA)
MetaData$depth <- as.factor(MetaData$depth)
MetaData$season <- as.factor(MetaData$month)

# sample names in both objects
samples_meta  <- MetaData$id
samples_count <- colnames(countData)
# find common samples
common_samples <- intersect(samples_meta, samples_count)
# subset and reorder count matrix
countData<- countData[, common_samples]
# reorder metadata to match countData
MetaData <- MetaData[match(common_samples, MetaData$id), ]
# must be TRUE
all(colnames(countData) == MetaData$id)

# filtering for low counts
keep <- rowSums(countData >= 10) >= 5

countData_filt <- countData[keep, ]

cat("Genes retained:", nrow(countData_filt), "\n")

# vst transformation

dds <- DESeqDataSetFromMatrix(
  countData = countData_filt,
  colData = MetaData,
  design = ~ site  + rRNA + RIN2 + condition 
)

dds <- estimateSizeFactors(dds)

vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# remove batch effect
mm <- model.matrix(~condition + site  + rRNA + RIN2 , colData(vsd))
treatment.design <- mm[,1:4]
batch.design <- mm[,-(1:4)]
mat <- limma::removeBatchEffect(vsd_mat, covariates=batch.design, design=treatment.design)
vsd_filtered <- mat

# PCA 
pca_input <- t(as.matrix(vsd_filtered))
pca <- prcomp(pca_input, scale. = FALSE) 

percent_var <- pca$sdev^2 / sum(pca$sdev^2)
percent_var_percent <- percent_var * 100

# Scree table
pca_variance_table <- tibble(
  PC = paste0("PC", seq_along(percent_var)),
  variance_explained = percent_var,
  percent_variance = percent_var_percent,
  cumulative_percent = cumsum(percent_var_percent)
)

head(pca_variance_table, 10)
#n_pcs <- which(cumsum(percent_var_percent) >= 70)[1]
n_pcs <- 2
pc_cols <- paste0("PC", 1:n_pcs)
pc_weights <- percent_var[1:n_pcs]

# Normalize weights 
pc_weights <- pc_weights / sum(pc_weights)

pc_weights


pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("id") %>%
  left_join(MetaData, by = "id")

# Keep a clean group column
pca_df <- pca_df %>%
  mutate(
    group = condition
  )

head(pca_df)

#### PERMANOVA ####
# PERMANOVA on the weighted distances

# PERMANOVA on the complete batch-corrected transcriptome profile
# Transpose 
global_dist <- dist(t(vsd_filtered)) 

# Run 
permanova_res <- adonis2(global_dist ~ depth * season, 
                         data = MetaData, 
                         permutations = 999)
print(permanova_res)
# Marginal testing 
permanova_terms <- adonis2(global_dist ~ depth * season, 
                           data = MetaData, 
                           permutations = 999,
                           by = "terms") 
print(permanova_terms)

pairwise_res <- pairwise.adonis2(global_dist ~ condition, data = MetaData, permutations = 999)

# Print pairwise comparisons
print(pairwise_res)

# Compute multivariate dispersion 
disp_res <- betadisper(global_dist, group = MetaData$condition)

# Test 
anova(disp_res)

permutest(disp_res, permutations = 999)
# betadisper not significant

# permanova for plot text
p_orig <- permanova_terms["depth", "Pr(>F)"]
p_trans <- permanova_terms["season", "Pr(>F)"]
p_text <- paste0("P_Depth = ", p_orig,"\nP_Season = ", p_trans) 

#### PCA distance to centroids ####
centroids_2d <- pca_df %>%
  group_by(condition) %>%
  summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  )

arrow_df <- tibble(
  shift = c("S to SS", "D to DD"),
  x = c(
    centroids_2d$PC1[centroids_2d$condition == "S"],
    centroids_2d$PC1[centroids_2d$condition == "D"]
  ),
  y = c(
    centroids_2d$PC2[centroids_2d$condition == "S"],
    centroids_2d$PC2[centroids_2d$condition == "D"]
  ),
  xend = c(
    centroids_2d$PC1[centroids_2d$condition == "SS"],
    centroids_2d$PC1[centroids_2d$condition == "DD"]
  ),
  yend = c(
    centroids_2d$PC2[centroids_2d$condition == "SS"],
    centroids_2d$PC2[centroids_2d$condition == "DD"]
  )
)

p_season_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 1.5) +
  #stat_ellipse(aes(group = condition), linewidth = 0.5, alpha = 0.35, show.legend = FALSE) +
  geom_point(
    data = centroids_2d,
    aes(x = PC1, y = PC2, color = condition),
    inherit.aes = FALSE,
    shape = 18,
    size = 2.3
  ) +
  geom_segment(
    data = arrow_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.1, "cm")),
    linewidth = 0.4,
    color = "black"
  ) +
  geom_text_repel(
    data = centroids_2d,
    aes(x = PC1, y = PC2, label = condition, color = condition),
    inherit.aes = FALSE,
    size = pt_to_mm(7),
    show.legend = FALSE
  ) +
  theme_classic() +
  labs(
    x = paste0("PC1: ", round(percent_var_percent[1], 1), "% variance"),
    y = paste0("PC2: ", round(percent_var_percent[2], 1), "% variance")
    #title = "PCA of VST-normalized host gene expression"
  ) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )+
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = p_text,
    hjust = -0.05,
    vjust = 1.2,
    size = pt_to_mm(6)
  )

p_season_pca

#### Variance-weighted distance function ####
weighted_pc_distance <- function(sample_scores, centroid_scores, weights) {
  sqrt(sum(weights * (sample_scores - centroid_scores)^2))
}

# Calculate SS and DD native/control centroids

centroids <- pca_df %>%
  filter(condition %in% c("S", "D")) %>%
  group_by(condition) %>%
  summarise(
    across(all_of(pc_cols), mean, na.rm = TRUE),
    .groups = "drop"
  )

S_centroid <- centroids %>%
  filter(condition == "S") %>%
  dplyr::select(all_of(pc_cols)) %>%
  as.numeric()

D_centroid <- centroids %>%
  filter(condition == "D") %>%
  dplyr::select(all_of(pc_cols)) %>%
  as.numeric()

# Calculate distances and indices

distance_season_df <- pca_df %>%
  rowwise() %>%
  mutate(
    dist_to_S = weighted_pc_distance(
      c_across(all_of(pc_cols)),
      S_centroid,
      pc_weights
    ),
    dist_to_D = weighted_pc_distance(
      c_across(all_of(pc_cols)),
      D_centroid,
      pc_weights
    )
  ) %>%
  ungroup() %>%
  mutate(
    winter_origin_centroid = case_when(
      depth == "S" ~ "S",
      depth == "D" ~ "D",
      TRUE ~ NA_character_
    ),
    
    distance_from_winter_origin = case_when(
      winter_origin_centroid == "S" ~ dist_to_S,
      winter_origin_centroid == "D" ~ dist_to_D,
      TRUE ~ NA_real_
    ),
    
    deep_like_index = dist_to_S / (dist_to_S + dist_to_D)
  )

# seasonal plasticity
pca_season_df <- distance_season_df %>%
  mutate(
    condition = factor(condition, levels = c("S", "SS", "D", "DD")),
    
    # seasonal displacement:
    # S and SS are measured from winter S centroid
    # D and DD are measured from winter D centroid
    Plasticity = distance_from_winter_origin,
    Seasonal_shift = distance_from_winter_origin
  )

# Check summary
pca_season_df %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean = mean(Plasticity, na.rm = TRUE),
    median = median(Plasticity, na.rm = TRUE),
    sd = sd(Plasticity, na.rm = TRUE),
    .groups = "drop"
  )

# ANOVA
fit_aov_season <- aov(Plasticity ~ condition, data = pca_season_df)

summary(fit_aov_season)
TukeyHSD(fit_aov_season)

# Assumption checks
shapiro.test(residuals(fit_aov_season))
car::leveneTest(Plasticity ~ condition, data = pca_season_df)
# not normal

# kruskal
kruskal.test(Plasticity ~ condition, data = pca_season_df)

pairwise_res <- pairwise.wilcox.test(
  pca_season_df$Plasticity,
  pca_season_df$condition,
  p.adjust.method = "BH"
)

pairwise_res

# letters

# Convert pairwise p-value matrix to letters
pw_p <- pairwise_res$p.value

# multcompLetters needs named vector
pw_vec <- pw_p[!is.na(pw_p)]
names(pw_vec) <- apply(which(!is.na(pw_p), arr.ind = TRUE), 1, function(i) {
  paste(rownames(pw_p)[i[1]], colnames(pw_p)[i[2]], sep = "-")
})

cld <- multcompLetters(pw_vec)

letter_df <- data.frame(
  condition = names(cld$Letters),
  Letter = cld$Letters
)

# plot
plas_stats <- pca_season_df %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean_p = mean(Plasticity, na.rm = TRUE),
    sd_p = sd(Plasticity, na.rm = TRUE),
    se_p = sd_p / sqrt(n),
    ci_val = 1.96 * se_p,
    lower = mean_p - ci_val,
    upper = mean_p + ci_val,
    max_group = max(Plasticity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(letter_df, by = "condition")

# Global y-position for significance letters
plas_stats <- plas_stats %>%
  mutate(label_y = 155)

p_season_shift <- ggplot(
  pca_season_df,
  aes(x = condition, y = Plasticity, color = condition)
) +
  # Raw data points
  geom_jitter(
    width = 0.15,
    size = 1.5,
    #alpha = 0.65,
    # shape = 1,
    # stroke = 1
  ) +
  # Mean point
  geom_point(
    data = plas_stats,
    aes(x = condition, y = mean_p),
    color = "black",
    size = 1.5,
    inherit.aes = FALSE
  ) +
  
  # 95% CI
  geom_errorbar(
    data = plas_stats,
    aes(x = condition, ymin = lower, ymax = upper),
    color = "black",
    width = 0,
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  
  # Significance letters
  geom_text(
    data = plas_stats,
    aes(x = condition, y = y_label, label = Letter),
    color = "black",
    size = pt_to_mm(6),
    inherit.aes = FALSE
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  theme_classic() +
  labs(
    x = "Treatment",
    y = "Distance from winter centroid"
    # title = "Weighted PCA transcriptomic plasticity"
  ) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    #axis.text.x = element_text(face = "bold"),
    legend.position = "none"
  )


p_season_shift


#### deep index ####

p_deepness <- ggplot(distance_season_df, aes(x = condition, y = deep_like_index, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(
    x = "Treatment",
    y = "Deep-like transcriptomic index",
    title = "Transcriptomic similarity to shallow and deep controls"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

p_deepness

# A. Deep-like index across all four treatments
model_deepness <- lm(deep_like_index ~ condition, data = distance_season_df)
anova(model_deepness)
TukeyHSD(aov(model_deepness))

# Assumption checks
shapiro.test(residuals(model_deepness))


#### plots ####
treatment_colors <- c(
  "S" = "#b495ff",
  "SS" = "#e16cb2",
  "SD" = "#4d68a4",
  "D" = "#9aced3",
  "DD" = "#58b799",
  "DS" = "#f97642"
)


combined <- (p_season_pca + p_season_shift)/ (p_pca + p_plas) +
  plot_annotation(tag_levels = 'A')& # Note the '&' instead of '+'
  theme(plot.tag = element_text(face = "bold", size = 9))&
  theme(plot.title = element_blank(), axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))&
  scale_color_manual(values = treatment_colors,
                     name = "Treatment") 
combined

ggsave(
  "/home/gospozha/haifa/cayman/manuscript/revision/figs/combined_fig5.pdf",
  plot = combined,
  width = 14.4,
  height = 11,
  units = "cm",
  device = "pdf",
  useDingbats = FALSE
)
ggsave(
  "/home/gospozha/haifa/cayman/manuscript/revision/figs/combined_fig5.jpg",
  plot = combined,
  width = 14.4,
  height = 11,
  units = "cm",
  device = "jpg"
)

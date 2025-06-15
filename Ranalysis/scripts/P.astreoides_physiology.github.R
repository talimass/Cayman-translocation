library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(car)
library(cowplot)
library(gridExtra)
library(dplyr)
library(multcompView)
library(EnvStats)
library(FSA)
library(emmeans)
library(ARTool)
library(vegan)
library(plotly)
library(rstatix)
library(tidyverse)
library(MASS)
library(patchwork)

# working directory
setwd("~/haifa/cayman/P.astreoides_physiology/github3/")

##### physiology #####
# opening file, checking, replacing characters to vectors when needed
physio <- read.csv('Physiology_P.astreoides2.csv', stringsAsFactors = F)
str(physio)
#View(physio)
physio$Treatment = factor(physio$Treatment, levels = c("S.T0", "SS", "SD", "D.T0", "DD", "DS"))
physio$Coral = factor(physio$Coral, levels = c("MF", "CC"))
physio$Month = factor(physio$Month, levels = c("Nov", "July"))
physio <- na.omit(physio)
# adding the same theme to each plot
mytheme = theme_bw()+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 8), axis.title = element_text(size = 9))

# analysis
# checking the effect of location 
vars <- physio %>% dplyr::select("cell.cm", "Protein.conc.ug.ml", "chl.surf")
adonis_result <- adonis2(vars ~ Coral, data = physio, method = "euclidean", permutations = 999)
print(adonis_result)
# no effect on vars, we can pool the samples

#### Protein.cm ####
# checking for normality and homoscedasticity
shapiro.test(physio$Protein.conc.ug.ml)  # data is not normal (p.value < 0.05)
leveneTest(Protein.conc.ug.ml~Treatment, d=physio) # heteroscedasticity of variance (p.value < 0.05)

#  transforming the data
shapiro.test(sqrt(physio$Protein.conc.ug.ml)) # data is normal after sqrt transformation
leveneTest(sqrt(Protein.conc.ug.ml)~Treatment, d=physio) # variance is ok

ggplot(physio, aes(x = sqrt(Protein.conc.ug.ml))) +
  geom_density() +
  mytheme

# performing anova + Tukey post-hoc 

model <- lm(data = physio, sqrt(Protein.conc.ug.ml) ~ Treatment)
# checking summary and plots
summary(model)
#plot(model) # Q-Q plots - checking for heteroscedasticity of residuals - look very good
# data is ok for anova
anova(model)

# posthoc test - pairwise comparisons (Tukey test)
posthoc <- TukeyHSD(aov(model))
posthoc$Treatment # the table of pairwise comparisons
# assigning letters to groups ("compact letter display")
letters <- multcompLetters4(aov(model), posthoc) 
# creating df of letters and their positions to add them to the plot
letters.df <- data.frame(letters$Treatment$Letters)
colnames(letters.df)[1] <- "Letter" 
letters.df$Treatment <- rownames(letters.df) 
placement <- physio %>% 
  group_by(Treatment) %>%
  summarise(quantile(Protein.conc.ug.ml)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement)

# boxplot with letters
# groups with the same letter are not significantly different
# groups that are significantly different get different letters
protein.anova <- ggplot(physio, aes(y = (Protein.conc.ug.ml), x = Treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  labs(y= "Protein concentration (ug/ml)", title = "Coral host protein concentration")+
  mytheme +
  guides(fill = FALSE)+
  scale_x_discrete(labels = c("10", "10→10", "10→40", "40", "40→40", "40→10"))+
  geom_boxplot(outlier.size = 0.5, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.4)  + 
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -1.25, vjust = -0.8, fontface = "italic")

protein.anova
#ggsave("prot.conc.anova.jpg", protein.anova, width = 4, height = 4)

#### cell/cm ####
# normality, homoscedasticity
shapiro.test(physio$cell.cm) # not normal  
leveneTest(cell.cm~Treatment*Coral,d=physio) # variance is ok

#transformation
shapiro.test(sqrt(physio$cell.cm)) # normal 
leveneTest(sqrt(cell.cm)~Treatment*Coral,d=physio) # variance is ok

ggplot(physio, aes(x = sqrt(cell.cm))) +
  geom_density() +
  mytheme

# anova + Tukey post-hoc 
model <- lm(data = physio, sqrt(cell.cm) ~ Treatment)
summary(model)  
#plot(model) # Q-Q plot is ok
anova(model) # Coral is not significant factor 

# Tukey posthoc test
posthoc <- TukeyHSD(aov(model)) 
posthoc$Coral # shows that CC is the same as MF
posthoc$Treatment # shows pairwise comparisons taking into account Coral in the model
letters <- multcompLetters4(aov(model), posthoc)
letters.df <- data.frame(letters$Treatment$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Treatment <- rownames(letters.df) 
placement <- physio %>%
  group_by(Treatment) %>%
  summarise(quantile(cell.cm, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement)

# boxplot with letters
cellcm.anova <- ggplot(physio, aes(y = (cell.cm), x = Treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  labs(y= ~cells ~x10^5 ~cm^-2, title = "Symbiont cell count per surface area")+
  mytheme +
  guides(fill = FALSE)+
  scale_x_discrete(labels = c("10", "10→10", "10→40", "40", "40→40", "40→10"))+
  geom_boxplot(outlier.size = 0.5, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.4)  + 
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -1.25, vjust = -0.8, fontface = "italic")
cellcm.anova
#ggsave("cellcm.anova.jpg", p.cellcm, width = 4, height = 4)

#### chl per surface ####
shapiro.test(physio$chl.surf) # data is not normal (p.value < 0.5)
leveneTest(chl.surf~Treatment,d=physio)  # heteroscedasticity of variance (p.value < 0.5)

# transforming the data (log, log10, sqrt, ^(1/3), exp,....)
shapiro.test(sqrt(physio$chl.surf)) # did not help
leveneTest((chl.surf)~Treatment,d=physio)  

# boxcox transformation
bc <-boxcox(lm(physio$chl.surf ~ physio$Treatment), lambda = seq(-2, 2, 0.1))
best_lambda <- bc$x[which.max(bc$y)]
abline(v = best_lambda, col = "red", lty = 2)
print(best_lambda)
physio$chl.surf.t <- (physio$chl.surf^best_lambda - 1) / best_lambda
shapiro.test((physio$chl.surf.t)) # normal
leveneTest((chl.surf.t)~Treatment,d=physio) #heteroscedascisity

# kruscal
# Kruskal-Wallis test for nonparametric data
kruskal.test(chl.surf ~ Treatment, data = physio) #  not significant
# not signif.
# posthoc: Dunn's Test with Bonferroni correction for p-values
dunn.res <- dunnTest(chl.surf ~ Treatment,
                     data=physio,
                     method="bonferroni")
diff <- dunn.res$res$P.adj < 0.05
Names <- gsub(" ", "", dunn.res$res$Comparison)
names(diff) <- Names
# compact letter display
letters <- multcompLetters(diff)
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter" 
letters.df$Treatment <- rownames(letters.df) 
placement <- physio %>%
  group_by(Treatment) %>%
  summarise(quantile(chl.surf, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 
# boxplot with letters

chl.surf<-ggplot(physio, aes(y = chl.surf, x = Treatment)) +
  labs(title = "Chlorophyll per surface area", y= ~chlorophyll[a] ~x10^-7 ~ cm^-2)+
  geom_boxplot(outlier.shape = NA, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  mytheme +
  guides(fill = FALSE)+
  scale_x_discrete(labels = c("10", "10→10", "10→40", "40", "40→40", "40→10"))+
  geom_boxplot(outlier.size = 0.5, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.4)  + 
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -1.25, vjust = -0.8, fontface = "italic")
chl.surf
##### photophysiology #####

fire.data <- read.csv("FIRe_tranplantation_exsitu_Nov2022_July2023.github.csv")
fire.data$Treatment = factor(fire.data$Treatment, levels = c("S.T0", "SS", "SD", "D.T0", "DD", "DS"))
fire.data$Coral <- as.factor(fire.data$Coral)
fire.data <- na.omit(fire.data)

# reduce dataset 
fire <- fire.data %>%
  drop_na() %>%
  group_by(coral.number, Treatment, Coral) %>%
  summarise(
    `Fv.Fm` = mean(`Fv.Fm`, na.rm = TRUE),
    Sigma = mean(Sigma, na.rm = TRUE),
    `Pmax.e.s` = mean(`Pmax.e.s`, na.rm = TRUE),
    p = mean(p, na.rm = TRUE),
    .groups = "drop"  # optional: ungroups the result
  )

# checking the effect of location
vars <- fire %>% dplyr::select(Fv.Fm, Sigma, Pmax.e.s, p)
adonis_result <- adonis2(vars ~ Coral, data = fire, method = "euclidean", permutations = 999)
print(adonis_result)
# no effect on vars

### Functional absorption cross-section of PSII ###

# checking for normality and homoscedasticity
shapiro.test((fire$Sigma))  # data is normal
leveneTest(Sigma~Treatment, d=fire) # heteroscedasticity of variance (p.value < 0.5)

# welch anova
oneway.test(data = fire, Sigma ~ Treatment, var.equal = TRUE)

# posthoc
games_results <- games_howell_test(fire, Sigma ~ Treatment, conf.level = 0.95, detailed = FALSE)

games_results$comparison <- paste0(games_results$group1,'-',games_results$group2)
# preparing pairwise p-values matrix
p_matrix <- games_results$p.adj.signif<0.05
names(p_matrix) <- games_results$comparison
# assigning letters to groups ("compact letter display")
letters <- multcompLetters(p_matrix)
# creating df of letters and their positions to add them to the plot
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Treatment <- rownames(letters.df) 
placement <- fire %>% 
  group_by(Treatment) %>%
  summarise(quantile(Sigma, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 

# boxplot with letters
# groups with the same letter are the same
# groups that are significantly different get different letters
sigma = ggplot(fire, aes(y = Sigma, x = Treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  labs(title="Functional absorption cross-section of PSII", y= "σPSII’(A2)") +
  mytheme +
  guides(fill = FALSE)+
  scale_x_discrete(labels = c("10", "10→10", "10→40", "40", "40→40", "40→10"))+
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

sigma

#### Quantum yield of photochemistry in PSII ####

# checking for normality and homoscedasticity

shapiro.test((fire$Fv.Fm))  # data is not normal (p.value < 0.05)
leveneTest(Fv.Fm~Treatment, d=fire) # heteroscedasticity of variance (p.value < 0.05)

# Kruskal-Wallis test for nonparametric data
kruskal.test(Fv.Fm ~ Treatment, data = fire) #  significant

# posthoc: Dunn's Test with Bonferroni correction for p-values
dunn.res <- dunnTest(Fv.Fm ~ Treatment,
                     data=fire,
                     method="bonferroni")

diff <- dunn.res$res$P.adj < 0.05
Names <- gsub(" ", "", dunn.res$res$Comparison)
names(diff) <- Names 
# assigning letters to groups ("compact letter display")
letters <- multcompLetters(diff)
# creating df of letters and their positions to add them to the plot
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Treatment <- rownames(letters.df) 
placement <- fire %>% 
  group_by(Treatment) %>%
  summarise(quantile(Fv.Fm, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 

# boxplot with letters
# groups with the same letter are the same
# groups that are significantly different get different letters
FvFm = ggplot(fire, aes(y = Fv.Fm, x = Treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  labs(title="Algal quantum yield of photochemistry in PSII", y="Fv’/Fm’") +
  mytheme +
  guides(fill = FALSE)+
  scale_x_discrete(labels = c("10", "10→10", "10→40", "40", "40→40", "40→10"))+
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

FvFm

#### Maximum photosynthetic rate ####
# checking for normality and homoscedasticity
shapiro.test((fire$Pmax.e.s))  # data is not normal (p.value < 0.05)
leveneTest(Pmax.e.s~Treatment, d=fire) # heteroscedasticity of variance (p.value < 0.05)

kruskal.test(Pmax.e.s ~ Treatment, data = fire) #  significant

# posthoc: Dunn's Test with Bonferroni correction for p-values
dunn.res <- dunnTest(Pmax.e.s ~ Treatment,
                     data=fire,
                     method="bonferroni")

diff <- dunn.res$res$P.adj < 0.05
Names <- gsub(" ", "", dunn.res$res$Comparison)
names(diff) <- Names 
# assigning letters to groups ("compact letter display")
letters <- multcompLetters(diff)
# creating df of letters and their positions to add them to the plot
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Treatment <- rownames(letters.df) 
placement <- fire %>% 
  group_by(Treatment) %>%
  summarise(quantile(Pmax.e.s, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 

# boxplot with letters
# groups with the same letter are the same
# groups that are significantly different get different letters
Pmax = ggplot(fire, aes(y = Pmax.e.s, x = Treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  labs(title="Maximum photosynthetic rate", y="Pmax (electron s-1 PSII-1)") +
  mytheme +
  guides(fill = FALSE)+
  scale_y_continuous(limits = c(60,130))+ 
  scale_x_discrete(labels = c("10", "10→10", "10→40", "40", "40→40", "40→10"))+
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

Pmax

#### connectivity parameter ####
# checking for normality and homoscedasticity
shapiro.test((fire$p))  # data is normal (p.value > 0.05)
leveneTest(p~Treatment, d=fire) # heteroscedasticity

# welch anova
oneway.test(data = fire, p ~ Treatment, var.equal = TRUE)

# posthoc
games_results <- games_howell_test(fire, p ~ Treatment, conf.level = 0.95, detailed = FALSE)

games_results$comparison <- paste0(games_results$group1,'-',games_results$group2)
# preparing pairwise p-values matrix
p_matrix <- games_results$p.adj.signif<0.05
names(p_matrix) <- games_results$comparison
# assigning letters to groups ("compact letter display")
letters <- multcompLetters(p_matrix)
# creating df of letters and their positions to add them to the plot
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Treatment <- rownames(letters.df) 
placement <- fire %>% 
  group_by(Treatment) %>%
  summarise(quantile(p, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 

# boxplot with letters
# groups with the same letter are the same
# groups that are significantly different get different letters
p = ggplot(fire, aes(y = p, x = Treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  labs(title="Connectivity parameter", y= "p") +
  mytheme +
  guides(fill = FALSE)+
  scale_x_discrete(labels = c("10", "10→10", "10→40", "40", "40→40", "40→10"))+
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

p

#### combining plots ####

# physiology
combined_phys <- (protein.anova | cellcm.anova | chl.surf) + plot_layout(heights = c(1, 1.2))

# Display it
combined_phys
ggsave("combined_phys.jpg", combined_phys, width = 12, height = 5)


# photochemistry
combined_phot <- (FvFm | sigma ) / (Pmax | p) + plot_layout(heights = c(1, 1.2))

# Display it
combined_phot
ggsave("combined_phot.jpg", combined_phot, width = 10, height = 10)


#### NMDS ####
joined_df <- full_join(physio, fire, by = c("coral.number", "Treatment", "Coral"))

joined_df  %>%
  drop_na() %>%
  dplyr::select(Treatment, Protein.conc.ug.ml, cell.cm, chl.surf, Fv.Fm, Sigma, Pmax.e.s, p)  -> nmds.physio
str(nmds.physio)

# Remove non-numeric columns (e.g., Sample, Treatment)
param_data <- nmds.physio[, -which(colnames(nmds.physio) %in% c("Treatment"))]

# Compute Bray-Curtis distance
bray_dist <- vegdist(param_data, method = "bray")

# Check distance matrix
as.matrix(bray_dist)
# Run NMDS with 2 dimensions (k=2)
nmds_result <- metaMDS(bray_dist, k = 2, trymax = 100)

# Check NMDS stress value (should be < 0.2 for a good fit)
nmds_result$stress
stress_value <- round(nmds_result$stress, 3) 
# Convert NMDS results to a dataframe
nmds_df <- as.data.frame(nmds_result$points)

# Add Treatment information
nmds_df$Treatment <- nmds.physio$Treatment
nmds_df$Treatment <- nmds.physio$Treatment
# Check final NMDS data
head(nmds_df)

# Fit environmental vectors (parameter influence)
env_fit <- envfit(nmds_result, param_data, perm = 999)

# Extract arrow coordinates
vectors <- as.data.frame(scores(env_fit, "vectors"))
vectors$Parameter <- rownames(vectors)  # Add parameter names
vectors$R2 <- round(env_fit$vectors$r, 2)  

nmdsp <-ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = Treatment)) +
  geom_point(size = 1) +                        # Sample points
  stat_ellipse(level = 0.95) +                  # Confidence ellipses
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +               # Vectors (arrows)
  geom_text(data = vectors, aes(x = NMDS1, y = NMDS2, label = Parameter), 
            vjust = -0.5, hjust = 0.5, color = "black", size = 3) +  
  annotate("text", x = min(nmds_result$points[,1]), y = min(nmds_result$points[,2]), 
           label = paste("Stress =", stress_value), size = 3, hjust = -2, color = "red") +
  scale_color_hue(labels = c("10",
                             "10→10","10→40",
                             "40",
                             "40→40", "40→10"))+
  theme_minimal() +
  labs(title = "NMDS of physiological parameters", x = "NMDS1", y = "NMDS2")
nmdsp
ggsave("nmds.all.jpg", nmdsp, width = 6, height = 6)



#### survival ####
surv <- read.csv('../NSF_transplant_survival.csv', stringsAsFactors = F)
str(surv)
surv$Treatment = factor(surv$Treatment, levels = c("SS", "SD", "DD", "DS"))
surv$Site = factor(surv$Site)
surv$Status = factor(surv$Status)
View(surv)
surv$Dead <- as.factor(ifelse(surv$Status == "D", 1, 0))

# For MF site
mf <- subset(surv, Site == "Marthas Finyard")

# Create contingency table: Dead vs Treatment
table_mf <- table(mf$Treatment, mf$Dead)

# Chi-squared test
chisq_mf <- chisq.test(table_mf)
chisq_mf

# For CC site
cc <- subset(surv, Site == "Coral City")

# Create contingency table: Dead vs Treatment
table_cc <- table(cc$Treatment, cc$Dead)

# Chi-squared test
chisq_cc <- chisq.test(table_cc)
chisq_cc

# Model for MF
model_mf <- glm(Dead ~ Treatment, data = mf, family = binomial)
emm_mf <- emmeans(model_mf, ~ Treatment)
pairs(emm_mf, adjust = "tukey")  # Post hoc comparison

# Model for CC
model_cc <- glm(Dead ~ Treatment, data = cc, family = binomial)
emm_cc <- emmeans(model_cc, ~ Treatment)
pairs(emm_cc, adjust = "tukey")

#  Prepare data with percentages
plot_data <- surv %>%
  group_by(Site, Treatment, Dead) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Site, Treatment) %>%
  mutate(Percent = 100 * N / sum(N))%>%
  mutate(Status = ifelse(Dead == 1, "Dead", "Alive"))

# Add significance stars on top of DD and DS bars
stars <- plot_data %>%
  filter(Treatment %in% c("DD", "DS")) %>%
  filter(Site %in% c("Marthas Finyard")) %>%
  group_by(Site, Treatment) %>%
  summarise(ypos = sum(Percent), .groups = "drop") %>%
  mutate(label = "*", ypos = ypos + 5)

# Step 3: Basic stacked barplot
s <- ggplot(plot_data, aes(x = Treatment, y = Percent, fill = Status)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site) +
  scale_fill_manual(values = c("Alive" = "aquamarine4", "Dead" = "coral")) +
  ylab("Percent") +
  theme_minimal() + 
  geom_text(aes(label = paste0(round(Percent, 1), "%")),
                              position = position_stack(vjust = 0.5),
                              size = 3, color = "white")+ 
  geom_text(data = stars, aes(x = Treatment, y = ypos, label = label),
                                                                    inherit.aes = FALSE, size = 6, fontface = "bold")
s
ggsave("survival.jpg", s, width = 6, height = 6)

save.image(file='physio4.RData')
load('physio4.RData')
#### growth ####

growth <- read.csv('../Alizarin-mark/growth.csv', stringsAsFactors = F)
str(growth)
growth$condition = factor(growth$condition, levels = c("SS", "SD", "DD", "DS"))
View(growth)

growthplot <- ggplot(growth, aes(y = size, x = condition)) +
  geom_boxplot(outlier.shape = NA, aes(fill = condition), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  labs(y= 'growth mm')+
  geom_jitter(position = position_jitter(width = .25), size = 0.5) +
  mytheme +
  guides(fill = FALSE)+
  scale_y_continuous(limits = c(1.3,3.3))+ 
  scale_x_discrete(labels = c("10→10", "10→40", "40→40", "40→10" ))
growthplot

ggsave("growth3.jpg", growthplot,  width = 6, height = 6)

shapiro.test(sqrt(growth$size)) # data is normal (p.value > 0.05)
leveneTest(sqrt(size)~condition,d=growth)  # no heteroscedasticity of variance (p.value > 0.05)

# performing anova + Tukey post-hoc 
growth <- na.omit(growth)
model <- lm(data = growth, sqrt(size) ~ condition)
# checking summary and plots
summary(model)

#plot(model) # Q-Q plots - checking for heteroscedasticity of residuals - look very good
# data is ok for anova
anova(model)

# posthoc test - pairwise comparisons (Tukey test)
posthoc <- TukeyHSD(aov(model))
posthoc$condition # the table of pairwise comparisons
# assigning letters to groups ("compact letter display")
letters <- multcompLetters4(aov(model), posthoc) 
# creating df of letters and their positions to add them to the plot
letters.df <- data.frame(letters$condition$Letters)
colnames(letters.df)[1] <- "Letter" 
letters.df$condition <- rownames(letters.df) 
placement <- growth %>% 
  group_by(condition) %>%
  summarise(quantile(size)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement)

# boxplot with letters
# groups with the same letter are not significantly different
# groups that are significantly different get different letters
growthplot <- ggplot(growth, aes(y = (size), x = condition)) +
  geom_boxplot(outlier.shape = NA, aes(fill = condition), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  labs(y= "mm", title = "Skeletal growth")+
  mytheme +
  geom_jitter(position = position_jitter(width = .25), size = 0.5) +
  guides(fill = FALSE)+
  scale_x_discrete(labels = c("10→10", "10→40", "40→40", "40→10" ))+
  geom_boxplot(outlier.size = 0.5, aes(fill = condition), fatten = 0.5, alpha = 0.7, lwd = 0.4)  + 
  geom_text(data = letters.df, aes(x = condition, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -1.25, vjust = -0.8, fontface = "italic")

growthplot
ggsave("growth3.jpg", growthplot,  width = 6, height = 6)


#### NMDS with growth (only for July) ####
names(growth) <- c("size", "coral.number", "Treatment", "Coral")
physio.growth <- left_join(joined_df, growth)

physio.growth  %>%
  drop_na() %>%
  dplyr::select(Treatment, Protein.conc.ug.ml, cell.cm, chl.surf, Fv.Fm, Sigma, Pmax.e.s, p, size)  -> nmds.physio
str(nmds.physio)

# Remove non-numeric columns (e.g., Sample, Treatment)
param_data <- nmds.physio[, -which(colnames(nmds.physio) %in% c("Treatment"))]

# Compute Bray-Curtis distance
bray_dist <- vegdist(param_data, method = "bray")

# Check distance matrix
as.matrix(bray_dist)
# Run NMDS with 2 dimensions (k=2)
nmds_result <- metaMDS(bray_dist, k = 2, trymax = 100)

# Check NMDS stress value (should be < 0.2 for a good fit)
nmds_result$stress
stress_value <- round(nmds_result$stress, 3) 
# Convert NMDS results to a dataframe
nmds_df <- as.data.frame(nmds_result$points)

# Add Treatment information
nmds_df$Treatment <- nmds.physio$Treatment
nmds_df$Treatment <- nmds.physio$Treatment
# Check final NMDS data
head(nmds_df)

# Fit environmental vectors (parameter influence)
env_fit <- envfit(nmds_result, param_data, perm = 999)

# Extract arrow coordinates
vectors <- as.data.frame(scores(env_fit, "vectors"))
vectors$Parameter <- rownames(vectors)  # Add parameter names
vectors$R2 <- round(env_fit$vectors$r, 2)  

nmdsp <-ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = Treatment)) +
  geom_point(size = 1) +                        # Sample points
  stat_ellipse(level = 0.95) +                  # Confidence ellipses
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +               # Vectors (arrows)
  geom_text(data = vectors, aes(x = NMDS1, y = NMDS2, label = Parameter), 
            vjust = -0.5, hjust = 0.5, color = "black", size = 3) +  
  annotate("text", x = min(nmds_result$points[,1]), y = min(nmds_result$points[,2]), 
           label = paste("Stress =", stress_value), size = 3, hjust = -2, color = "red") +
  scale_color_hue(labels = c("10→10","10→40",
                             "40→40", "40→10"))+
  theme_minimal() +
  labs(title = "NMDS of physiological parameters (including growth)", x = "NMDS1", y = "NMDS2")
nmdsp
ggsave("nmds.july.growth.jpg", nmdsp, width = 6, height = 6)

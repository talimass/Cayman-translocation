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
library(pairwiseAdonis)
library(plotly)
library(rstatix)
library(tidyverse)
library(MASS)
library(patchwork)
library(goft)
library(boot)
library(ggpattern)
library(lme4)
library(lmerTest)
library(fitdistrplus)
# working directory
setwd("~/haifa/cayman/P.astreoides_physiology/github3/")
save.image("physio.revision.R")
load("physio.revision.R")
##### physiology #####
# opening file, checking, replacing characters to vectors when needed
physio <- read.csv('Physiology_P.astreoides2_revision.csv', stringsAsFactors = F)
str(physio)
#View(physio)
physio$Treatment = factor(physio$Treatment, levels = c("S", "SS", "SD", "D", "DD", "DS"))
physio$Coral = factor(physio$Coral, levels = c("MF", "CC"))
physio$Month = factor(physio$Month, levels = c("Nov", "July"))
physio$colony_id <- paste(physio$colony, physio$origin, physio$Coral, sep = "_")
physio$coral_id <- paste(physio$coral.number, physio$origin, physio$Coral, sep = "_")
physio <- na.omit(physio)
# adding the same theme to each plot
mytheme = theme_bw()+
  theme(axis.title.x = element_blank(),plot.title = element_text(size = 10),
        axis.text = element_text(colour = "black", size = 8), axis.title = element_text(size = 9))

# defining custom colors and patterns
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


# analysis
# checking the effect of location 
vars <- physio %>% dplyr::select("cell.cm", "Protein.conc.ug.ml", "chl.cell")
adonis_result <- adonis2(vars ~ Coral, data = physio, method = "euclidean", permutations = 999)
print(adonis_result)
# no effect on vars

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

# lmer
# mixed model
model <- lmer(sqrt(Protein.conc.ug.ml) ~ Treatment + (1 | colony_id), data = physio)

summary(model)
anova(model)
nobs(model)
# diagnostics
plot(model)
table(physio$Treatment)
levels(physio$Treatment)
qqnorm(resid(model))
qqline(resid(model))

# posthoc Tukey comparisons
emm <- emmeans(model, ~ Treatment)

posthoc <- pairs(emm, adjust = "tukey")
# posthoc test
posthoc <- emmeans(model, pairwise ~ Treatment, adjust = "tukey")
p_values <- as.data.frame(posthoc$contrasts)
d <- p_values$p.value < 0.05
Names <- gsub(" ", "", p_values$contrast)
names(d) <- Names
# compact letter display
letters <- multcompLetters(d)
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Treatment <- rownames(letters.df) 
placement <- physio %>% 
  group_by(Treatment) %>%
  summarise(quantile(Protein.conc.ug.ml, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 

# boxplot with letters
protein.anova <- ggplot(physio, aes(y = (Protein.conc.ug.ml), x = Treatment)) +
  geom_boxplot_pattern(
    aes(fill = Treatment, pattern = Treatment), outlier.size = 0.5, fatten = 0.5, alpha = 1,
    lwd = 0.4, pattern_fill = "gray", pattern_density = 0.08, pattern_spacing = 0.02) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = fill_patterns) +
  labs(y= "Protein concentration (ug/ml)", title = "Coral host protein concentration")+
  mytheme +
  guides(fill = FALSE, pattern = FALSE) +
  #scale_x_discrete(labels = c("10", "10→10", "10→40", "40", "40→40", "40→10")) +
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

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

# lmer
# mixed model
model <- lmer(sqrt(cell.cm) ~ Treatment + (1 | colony_id), data = physio)

summary(model)
anova(model)
nobs(model)
# diagnostics
plot(model)
qqnorm(resid(model))
qqline(resid(model))

# posthoc Tukey comparisons
emm <- emmeans(model, ~ Treatment)

posthoc <- pairs(emm, adjust = "tukey")
# posthoc test
posthoc <- emmeans(model, pairwise ~ Treatment, adjust = "tukey")
p_values <- as.data.frame(posthoc$contrasts)
d <- p_values$p.value < 0.05
Names <- gsub(" ", "", p_values$contrast)
names(d) <- Names
# compact letter display
letters <- multcompLetters(d)
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Treatment <- rownames(letters.df) 
placement <- physio %>% 
  group_by(Treatment) %>%
  summarise(quantile(cell.cm, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 

# boxplot with letters

cellcm.anova <- ggplot(physio, aes(y = (cell.cm), x = Treatment)) +
  geom_boxplot_pattern(
    aes(fill = Treatment, pattern = Treatment), outlier.size = 0.5, fatten = 0.5, alpha = 1,
    lwd = 0.4, pattern_fill = "gray", pattern_density = 0.08, pattern_spacing = 0.02) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = fill_patterns) +
  labs(y= ~cells ~x10^5 ~cm^-2, title = "Symbiont cell count per surface area")+
  mytheme +
  guides(fill = FALSE, pattern = FALSE) +
  # scale_x_discrete(labels = c("10", "10→10", "10→40", "40", "40→40", "40→10")) +
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

cellcm.anova


#### chl per cell ####

shapiro.test(physio$chl.cell) # data is not normal (p.value < 0.5)
leveneTest(chl.cell~Treatment,d=physio)  # heteroscedasticity of variance (p.value < 0.5)

# transforming the data (log, log10, sqrt, ^(1/3), exp,....)
shapiro.test(sqrt(physio$chl.cell)) # did not help

# checking if data follows gamma or inverse gaussian distribution
gamma_test(physio$chl.cell) # doesn't fit, p.value < 0.05
# inverse gaussian
ig_test(physio$chl.cell, method = "transf") # fit

# lmer

# GLMM with inverse gaussian family
# trying two factors
model <- glmer(
  chl.cell ~ Treatment + (1 | colony_id),
  data = physio,
  family = inverse.gaussian(link = "log")
)
#model1 <- glm(chl.cell ~ Treatment + Coral, data = physio, family = inverse.gaussian(link = "log"))
#AIC(model, model1) # model is better
##model2 <- glm(chl.cell ~ Treatment, data = physio, family = inverse.gaussian(link = "log"))
#AIC(model, model2) # model 2 is better
summary(model)
plot(model)
anova(model)
car::Anova(model, type = 3)
nobs(model)

# reduced/null model
model_reduced <- glmer(
  chl.cell ~ 1 + (1 | colony_id),
  data = physio,
  family = inverse.gaussian(link = "log")
)

# likelihood-ratio test
anova(model_reduced, model, test = "Chisq")


# posthoc test
posthoc <- emmeans(model, pairwise ~ Treatment, adjust = "tukey")
p_values <- as.data.frame(posthoc$contrasts)
d <- p_values$p.value < 0.05
Names <- gsub(" ", "", p_values$contrast)
names(d) <- Names
# compact letter display
letters <- multcompLetters(d)
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Treatment <- rownames(letters.df) 
placement <- physio %>% 
  group_by(Treatment) %>%
  summarise(quantile(chl.cell, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 
letters.df[6, 3] <- 21 # changing this position because of the outlier

# boxplot with letters

chlcell = ggplot(physio, aes(y = chl.cell, x = Treatment)) +
  geom_boxplot_pattern(
    aes(fill = Treatment, pattern = Treatment), outlier.size = 0.5, fatten = 0.5, alpha = 1,
    lwd = 0.4, pattern_fill = "gray", pattern_density = 0.08, pattern_spacing = 0.02) +
  scale_y_continuous(limits = c(0,53))+ 
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = fill_patterns) +
  labs(title = "Chlorophyll per algae cell", y= ~chlorophyll[a] ~x10^-7 ~cell^-1)+
  mytheme +
  guides(fill = FALSE, pattern = FALSE) +
  #scale_x_discrete(labels = c("10", "10→10", "10→40", "40", "40→40", "40→10")) +
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.7, fontface = "italic")

chlcell
##### photophysiology #####

fire.data <- read.csv("FIRe_tranplantation_exsitu_Nov2022_July2023.github_revision.csv")
fire.data$Treatment = factor(fire.data$Treatment, levels = c("S", "SS", "SD", "D", "DD", "DS"))
fire.data$Coral <- as.factor(fire.data$Coral)
fire.data$coral_id <- paste(fire.data$coral.number, fire.data$origin, fire.data$Coral, sep = "_")
fire.data <- na.omit(fire.data)

# reduce dataset 
fire <- fire.data %>%
  drop_na() %>%
  group_by(coral.number, Treatment, Coral, coral_id) %>%
  summarise(
    `Fv.Fm` = mean(`Fv.Fm`, na.rm = TRUE),
    Sigma = mean(Sigma, na.rm = TRUE),
    `Pmax.e.s` = mean(`Pmax.e.s`, na.rm = TRUE),
    p = mean(p, na.rm = TRUE),
    .groups = "drop"  # optional: ungroups the result
  )

fire <- left_join(
  fire,
  physio %>%
    dplyr::select(coral_id, colony_id) %>%
    distinct(),
  by = "coral_id"
)

fire[65, 9] <- "5_D_MF"
fire[72, 9] <- "5_S_MF"
fire[73, 9] <- "1_S_CC"
fire[89, 9] <-  "3_D_CC"
fire[90, 9] <- "5_D_CC"
  
# checking the effect of location
vars <- fire %>% dplyr::select(Fv.Fm, Sigma, Pmax.e.s, p)
adonis_result <- adonis2(vars ~ Coral, data = fire, method = "euclidean", permutations = 999)
print(adonis_result)
# no effect on vars

### Functional absorption cross-section of PSII ###

# checking for normality and homoscedasticity
shapiro.test((fire$Sigma))  # data is normal
leveneTest(Sigma~Treatment, d=fire) # heteroscedasticity of variance (p.value < 0.5)


# lmer
library(nlme)
model_sigma <- lme(
  Sigma ~ Treatment,
  random = ~1 | colony_id,
  weights = varIdent(form = ~1 | Treatment),
  data = fire,
  method = "REML"
)

summary(model_sigma)
anova(model_sigma)

# diagnostics
plot(model_sigma)
qqnorm(resid(model_sigma, type = "normalized"))
qqline(resid(model_sigma, type = "normalized"))

# posthoc Tukey comparisons


posthoc <- emmeans(model_sigma, pairwise ~ Treatment, adjust = "tukey")
p_values <- as.data.frame(posthoc$contrasts)
d <- p_values$p.value < 0.05
Names <- gsub(" ", "", p_values$contrast)
names(d) <- Names
# compact letter display
letters <- multcompLetters(d)
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Treatment <- rownames(letters.df) 
placement <- fire %>% 
  group_by(Treatment) %>%
  summarise(quantile(Sigma, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 

sigma = ggplot(fire, aes(y = Sigma, x = Treatment)) +
  geom_boxplot_pattern(
    aes(fill = Treatment, pattern = Treatment), outlier.size = 0.5, fatten = 0.5, alpha = 1,
    lwd = 0.4, pattern_fill = "gray", pattern_density = 0.08, pattern_spacing = 0.02) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = fill_patterns) +
  labs(title="Functional absorption cross-section of PSII", y= "σPSII’(A2)") +
  mytheme +
  guides(fill = FALSE, pattern = FALSE) +
  #scale_x_discrete(labels = c("10", "10→10", "10→40", "40", "40→40", "40→10")) +
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

sigma

#### Quantum yield of photochemistry in PSII ####

# checking for normality and homoscedasticity

shapiro.test((fire$Fv.Fm))  # data is not normal (p.value < 0.05)
leveneTest(Fv.Fm~Treatment, d=fire) # heteroscedasticity of variance (p.value < 0.05)



# lmer
model_fvfm <- lme(
  Fv.Fm ~ Treatment,
  random = ~1 | colony_id,
  weights = varIdent(form = ~1 | Treatment),
  data = fire,
  method = "REML"
)

qqnorm(resid(model_fvfm, type = "normalized"))
qqline(resid(model_fvfm, type = "normalized"))

hist(resid(model_fvfm, type = "normalized"))
plot(model_fvfm)
# posthoc Tukey comparisons
anova(model_fvfm)

posthoc <- emmeans(model_fvfm, pairwise ~ Treatment, adjust = "tukey")
p_values <- as.data.frame(posthoc$contrasts)
d <- p_values$p.value < 0.05
Names <- gsub(" ", "", p_values$contrast)
names(d) <- Names
# compact letter display
letters <- multcompLetters(d)
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
  geom_boxplot_pattern(
    aes(fill = Treatment, pattern = Treatment), outlier.size = 0.5, fatten = 0.5, alpha = 1,
    lwd = 0.4, pattern_fill = "gray", pattern_density = 0.08, pattern_spacing = 0.02) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = fill_patterns) +
  labs(title="Algal quantum yield of photochemistry in PSII", y="Fv’/Fm’") +
  mytheme +
  guides(fill = FALSE, pattern = FALSE) +
  #scale_x_discrete(labels = c("10", "10→10", "10→40", "40", "40→40", "40→10")) +
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

FvFm

#### Maximum photosynthetic rate ####
# checking for normality and homoscedasticity
shapiro.test((fire$Pmax.e.s))  # data is not normal (p.value < 0.05)
leveneTest(Pmax.e.s~Treatment, d=fire) # heteroscedasticity of variance (p.value < 0.05)

#lmer


descdist(fire$Pmax.e.s, boot = 1000)

fire %>%
  arrange(desc(Pmax.e.s)) %>%
  dplyr::select(coral_id, colony_id, Treatment, Pmax.e.s) %>%
  head(10)

model_pmax <- lme(
  log(Pmax.e.s) ~ Treatment,
  random = ~1 | colony_id,
  weights = varIdent(form = ~1 | Treatment),
  data = fire,
  method = "REML"
)
qqnorm(resid(model_pmax, type = "normalized"))
qqline(resid(model_pmax, type = "normalized"))

plot(model_pmax)

anova(model_pmax)

posthoc <- emmeans(model_pmax, pairwise ~ Treatment, adjust = "tukey")
p_values <- as.data.frame(posthoc$contrasts)
d <- p_values$p.value < 0.05
Names <- gsub(" ", "", p_values$contrast)
names(d) <- Names
# compact letter display
letters <- multcompLetters(d)
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Treatment <- rownames(letters.df) 
placement <- fire %>% 
  group_by(Treatment) %>%
  summarise(quantile(Pmax.e.s, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 


Pmax = ggplot(fire, aes(y = Pmax.e.s, x = Treatment)) +
  geom_boxplot_pattern(
    aes(fill = Treatment, pattern = Treatment), outlier.size = 0.5, fatten = 0.5, alpha = 1,
    lwd = 0.4, pattern_fill = "gray", pattern_density = 0.08, pattern_spacing = 0.02) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = fill_patterns) +
  ylim(60, 130) +
  labs(title="Maximum photosynthetic rate", y="Pmax (electron s-1 PSII-1)") +
  mytheme +
  guides(fill = FALSE, pattern = FALSE) +
  #scale_x_discrete(labels = c("10", "10→10", "10→40", "40", "40→40", "40→10")) +
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

Pmax
#### connectivity parameter ####
# checking for normality and homoscedasticity
shapiro.test((fire$p))  # data is normal (p.value > 0.05)
leveneTest(p~Treatment, d=fire) # heteroscedasticity


# lmer
model_fvfm <- lme(
  p ~ Treatment,
  random = ~1 | colony_id,
  weights = varIdent(form = ~1 | Treatment),
  data = fire,
  method = "REML"
)

anova(model_fvfm)
qqnorm(resid(model_fvfm, type = "normalized"))
qqline(resid(model_fvfm, type = "normalized"))

hist(resid(model_fvfm, type = "normalized"))
plot(model_fvfm)
# posthoc Tukey comparisons

posthoc <- emmeans(model_fvfm, pairwise ~ Treatment, adjust = "tukey")
p_values <- as.data.frame(posthoc$contrasts)
d <- p_values$p.value < 0.05
Names <- gsub(" ", "", p_values$contrast)
names(d) <- Names
# compact letter display
letters <- multcompLetters(d)
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Treatment <- rownames(letters.df) 
placement <- fire %>% 
  group_by(Treatment) %>%
  summarise(quantile(p, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 

p = ggplot(fire, aes(y = p, x = Treatment)) +
  geom_boxplot_pattern(
    aes(fill = Treatment, pattern = Treatment), outlier.size = 0.5, fatten = 0.5, alpha = 1,
    lwd = 0.4, pattern_fill = "gray", pattern_density = 0.08, pattern_spacing = 0.02) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = fill_patterns) +
  labs(title="Connectivity parameter", y= "p") +
  mytheme +
  guides(fill = FALSE, pattern = FALSE) +
  #scale_x_discrete(labels = c("10", "10→10", "10→40", "40", "40→40", "40→10")) +
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

p


#### corallite density ####
corallite <- read.csv('../Polyp_size.csv', stringsAsFactors = F)
str(corallite)

# Set treatment order
corallite$Treatment <- factor(corallite$Treatment, levels = c("SS", "SD", "DD", "DS"))

# Make coral number a factor
corallite$Coral.number <- factor(corallite$Coral.number)

corallite_dens <- na.omit(corallite)

m_pois <- glm(
  Polyp.density.4_95.2mm ~ Treatment,
  data = corallite_dens,
  family = poisson(link = "log")
)

summary(m_pois)
nobs(m_pois)
dispersion <- sum(residuals(m_pois, type = "pearson")^2) / df.residual(m_pois)
dispersion


m_null <- glm(
  Polyp.density.4_95.2mm ~ 1,
  data = corallite_dens,
  family = poisson(link = "log")
)

anova(m_null, m_pois, test = "Chisq")

emm_pois <- emmeans(m_pois, ~ Treatment, type = "response")
emm_pois
pairs(emm_pois)

# boxplot with letters
# groups with the same letter are the same
# groups that are significantly different get different letters

corallite_density = ggplot(corallite_dens, aes(y = Polyp.density.4_95.2mm, x = Treatment)) +
  geom_boxplot_pattern(
    aes(fill = Treatment, pattern = Treatment), outlier.size = 0.5, fatten = 0.5, alpha = 1,
    lwd = 0.4, pattern_fill = "gray", pattern_density = 0.08, pattern_spacing = 0.02) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = fill_patterns) +
  labs(y="Number of polyps per 5 mm²") +
  mytheme +
  guides(fill = FALSE, pattern = FALSE) 
corallite_density


#### corallite diameter ####

corallite <- read.csv('../Polyp_size_260526.csv', stringsAsFactors = F)
str(corallite)

# Set treatment order
corallite$Treatment <- factor(corallite$Treatment, levels = c("SS", "SD", "DD", "DS"))

# Make coral number a factor
corallite$Coral.number <- factor(corallite$Coral.number)

hist(corallite$Corallite.Diameter.mm)
shapiro.test(corallite$Corallite.Diameter.mm) # not normal but we'll check residuals after
leveneTest(Corallite.Diameter.mm~Treatment, d=corallite) # heteroscedasticity

shapiro.test(log(corallite$Corallite.Diameter.mm)) # not normal but we'll check residuals after
leveneTest(log(Corallite.Diameter.mm)~Treatment, d=corallite) # heteroscedasticity

# Linear Mixed-Effects Models, diff polyps = random effect

model_corallite_log <- lmer(
  log(Corallite.Diameter.mm) ~ Treatment + (1 | Coral.number),
  data = corallite
)

shapiro.test(resid(model_corallite_log))

qqnorm(resid(model_corallite_log))
qqline(resid(model_corallite_log))

plot(model_corallite_log)
anova(model_corallite_log)
nobs(model_corallite_log)
# not significant


corallite_plot <- ggplot(corallite, aes(y = (Corallite.Diameter.mm), x = Treatment)) +
  geom_boxplot_pattern(
    aes(fill = Treatment, pattern = Treatment), outlier.size = 0.5, fatten = 0.5, alpha = 1,
    lwd = 0.4, pattern_fill = "gray", pattern_density = 0.08, pattern_spacing = 0.02) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = fill_patterns) +
  #labs(y= "mm", title = "Skeletal extension rate")+
  labs(y= "mm")+
  mytheme +
  guides(fill = FALSE, pattern = FALSE) 
corallite_plot


# ggsave("5/corallite_diameter.jpg", corallite_plot_mean,  width = 6, height = 6)
# ggsave("5/corallite_diameter.pdf", corallite_plot_mean,  width = 6, height = 6)
# ggsave("5/corallite_diameter3.jpg", corallite_plot,  width = 6, height = 6)
# ggsave("5/corallite_diameter3.pdf", corallite_plot,  width = 6, height = 6)
# 


corallite_combined <-  (corallite_plot + corallite_density) +
  plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(face = "bold")) 

corallite_combined <-
  corallite_combined&
  theme(plot.title = element_blank(), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        text = element_text(size = 7))



# ggsave("5/corallite_all.jpg", corallite_combined,  width = 8, height = 4)
# ggsave("5/corallite_all.pdf", corallite_combined,  width = 8, height = 4)


ggsave(
  "/home/gospozha/haifa/cayman/manuscript/revision/figs/combined_figS5.pdf",
  plot = corallite_combined,
  width = 18.4,
  height = 9,
  units = "cm",
  device = "pdf",
  useDingbats = FALSE
)
ggsave(
  "/home/gospozha/haifa/cayman/manuscript/revision/figs/combined_figS5.jpg",
  plot = corallite_combined,
  width = 18.4,
  height = 9,
  units = "cm",
  device = "jpg"
)
#### growth ####

growth <- read.csv('../Alizarin-mark/growth_revision.csv', stringsAsFactors = F)
str(growth)
growth$condition = factor(growth$condition, levels = c("SS", "SD", "DD", "DS"))
#growth$Coral <- as.factor(growth$Coral)
growth$coral_id <- paste(growth$sample, growth$origin, growth$site, sep = "_")

growth <- left_join(
  growth,
  physio %>%
    dplyr::select(coral_id, colony_id) %>%
    distinct(),
  by = "coral_id"
)


shapiro.test(sqrt(growth$size)) # data is normal (p.value > 0.05)
leveneTest(sqrt(size)~condition,d=growth)  # no heteroscedasticity of variance (p.value > 0.05)

# lmer
# mixed model
model <- lmer(sqrt(size) ~ condition + (1 | colony_id), data = growth)

summary(model)
anova(model)
nobs(model)
# diagnostics
plot(model)
qqnorm(resid(model))
qqline(resid(model))

# posthoc Tukey comparisons
emm <- emmeans(model, ~ condition)

posthoc <- pairs(emm, adjust = "tukey")
# posthoc test
posthoc <- emmeans(model, pairwise ~ condition, adjust = "tukey")
p_values <- as.data.frame(posthoc$contrasts)
d <- p_values$p.value < 0.05
Names <- gsub(" ", "", p_values$contrast)
names(d) <- Names
# compact letter display
letters <- multcompLetters(d)
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$condition <- rownames(letters.df) 
placement <- growth %>% 
  group_by(condition) %>%
  summarise(quantile(size, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 

# boxplot with letters
# groups with the same letter are not significantly different
# groups that are significantly different get different letters

growthplot <- ggplot(growth, aes(y = (size), x = condition)) +
  geom_boxplot_pattern(
    aes(fill = condition, pattern = condition), outlier.size = 0.5, fatten = 0.5, alpha = 1,
    lwd = 0.4, pattern_fill = "gray", pattern_density = 0.08, pattern_spacing = 0.02) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = fill_patterns) +
  #labs(y= "mm", title = "Skeletal extension rate")+
  labs(y= "mm")+
  mytheme +
  guides(fill = FALSE, pattern = FALSE) +
  # scale_x_discrete(labels = c("10→10", "10→40", "40→40", "40→10" ))+
  geom_text(data = letters.df, aes(x = condition, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

growthplot


ggsave("5/growth.jpg", growthplot,  width = 6, height = 6)
ggsave("5/growth.pdf", growthplot,  width = 6, height = 6)

# growth in cm/year
year_fraction <- 230 / 365
year_fraction
growth$growth_cm_yr <- (growth$size / 10) / year_fraction

shapiro.test(sqrt(growth$growth_cm_yr)) # data is normal (p.value > 0.05)
leveneTest(sqrt(growth_cm_yr)~condition,d=growth)  # no heteroscedasticity of variance (p.value > 0.05)

# performing anova + Tukey post-hoc 
growth <- na.omit(growth)
model <- lm(data = growth, sqrt(growth_cm_yr) ~ condition)
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
  summarise(quantile(growth_cm_yr)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement)

# boxplot with letters
# groups with the same letter are not significantly different
# groups that are significantly different get different letters

growthplot <- ggplot(growth, aes(y = (growth_cm_yr), x = condition)) +
  geom_boxplot_pattern(
    aes(fill = condition, pattern = condition), outlier.size = 0.5, fatten = 0.5, alpha = 1,
    lwd = 0.4, pattern_fill = "gray", pattern_density = 0.08, pattern_spacing = 0.02) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = fill_patterns) +
  #labs(y= "mm", title = "Skeletal extension rate")+
  labs(y=  expression(cm~yr^{-1})) +
  mytheme +
  guides(fill = FALSE, pattern = FALSE) +
  # scale_x_discrete(labels = c("10→10", "10→40", "40→40", "40→10" ))+
  geom_text(data = letters.df, aes(x = condition, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

growthplot
ggsave("5/growth_cm_year.jpg", growthplot,  width = 4, height = 4)
ggsave("5/growth_cm_year.pdf", growthplot,  width = 4, height = 4)

# average and sd

summary_growth <- growth %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean_growth = mean(growth_cm_yr, na.rm = TRUE),
    se_growth = sd(growth_cm_yr, na.rm = TRUE) / sqrt(n)
  )

summary_growth


# band in mm, combine
model <- lmer(sqrt(band) ~ condition + (1 | colony_id), data = growth)

summary(model)
anova(model)
nobs(model)
# diagnostics
plot(model)
qqnorm(resid(model))
qqline(resid(model))

# posthoc Tukey comparisons
emm <- emmeans(model, ~ condition)

posthoc <- pairs(emm, adjust = "tukey")
# posthoc test
posthoc <- emmeans(model, pairwise ~ condition, adjust = "tukey")
p_values <- as.data.frame(posthoc$contrasts)
d <- p_values$p.value < 0.05
Names <- gsub(" ", "", p_values$contrast)
names(d) <- Names
# compact letter display
letters <- multcompLetters(d)
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$condition <- rownames(letters.df) 
placement <- growth %>% 
  group_by(condition) %>%
  summarise(quantile(size, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 
model <- lm(data = growth, sqrt(band) ~ condition)
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
  summarise(quantile(band)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement)

# boxplot with letters
# groups with the same letter are not significantly different
# groups that are significantly different get different letters

growthplot <- ggplot(growth, aes(y = (band), x = condition)) +
  geom_boxplot_pattern(
    aes(fill = condition, pattern = condition), outlier.size = 0.5, fatten = 0.5, alpha = 1,
    lwd = 0.4, pattern_fill = "gray", pattern_density = 0.08, pattern_spacing = 0.02) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = fill_patterns) +
  #labs(y= "mm", title = "Skeletal extension rate")+
  labs(y= "mm")+
  mytheme +
  guides(fill = FALSE, pattern = FALSE) +
  # scale_x_discrete(labels = c("10→10", "10→40", "40→40", "40→10" ))+
  geom_text(data = letters.df, aes(x = condition, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

growthplot

# band avg
summary_growth <- growth %>%
  summarise(
    n = n(),
    mean_growth = mean(band, na.rm = TRUE),
    se_growth = sd(band, na.rm = TRUE) / sqrt(n)
  )

summary_growth

# band in cm/year
year_fraction <- 230 / 365
year_fraction
growth$band_cm_yr <- (growth$rest / 10) / year_fraction

shapiro.test(sqrt(growth$band_cm_yr)) # data is normal (p.value > 0.05)
leveneTest(sqrt(band_cm_yr)~condition,d=growth)  # no heteroscedasticity of variance (p.value > 0.05)

# performing anova + Tukey post-hoc 
growth <- na.omit(growth)
model <- lm(data = growth, sqrt(band_cm_yr) ~ condition)
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
  summarise(quantile(band_cm_yr)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement)

# boxplot with letters
# groups with the same letter are not significantly different
# groups that are significantly different get different letters

growthplot <- ggplot(growth, aes(y = (band_cm_yr), x = condition)) +
  geom_boxplot_pattern(
    aes(fill = condition, pattern = condition), outlier.size = 0.5, fatten = 0.5, alpha = 1,
    lwd = 0.4, pattern_fill = "gray", pattern_density = 0.08, pattern_spacing = 0.02) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = fill_patterns) +
  #labs(y= "mm", title = "Alizarin")+
  labs(y=  expression(cm~yr^{-1})) +
  mytheme +
  guides(fill = FALSE, pattern = FALSE) +
  # scale_x_discrete(labels = c("10→10", "10→40", "40→40", "40→10" ))+
  geom_text(data = letters.df, aes(x = condition, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

growthplot
ggsave("5/band_cm_year.jpg", growthplot,  width = 4, height = 4)
ggsave("5/band_cm_year.pdf", growthplot,  width = 4, height = 4)



#### combining plots ####

combined_phys <- (
  protein.anova | cellcm.anova | chlcell) / (FvFm | Pmax )/( sigma | p)  +
  plot_annotation(tag_levels = "A")& 
  theme(plot.tag = element_text(face = "bold", size = 9)) 

combined_phys <- 
  combined_phys &
  theme(plot.title = element_blank(), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        text = element_text(size = 7))

combined_phys


ggsave(
  "/home/gospozha/haifa/cayman/manuscript/revision/figs/combined_figS6.pdf",
  plot = combined_phys,
  width = 18.4,
  height = 19.9,
  units = "cm",
  device = "pdf",
  useDingbats = FALSE
)
ggsave(
  "/home/gospozha/haifa/cayman/manuscript/revision/figs/combined_figS6.jpg",
  plot = combined_phys,
  width = 18.4,
  height = 19.9,
  units = "cm",
  device = "jpg"
)
# display it
# ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/combined3_notitles_revision.jpg", combined_phys, width = 9, height =11)
# ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/combined3_notitles_revision.pdf", combined_phys, width = 9, height =11)


#### normalizing data before pca ####
joined_df <- full_join(physio, fire, by = c("coral.number", "Treatment", "Coral"))

joined_df  %>%
  drop_na() %>%
  dplyr::select(Protein.conc.ug.ml, cell.cm, chl.cell, Fv.Fm, Sigma, Pmax.e.s, p)  %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  ) -> dat_long

ggplot(dat_long, aes(x = value)) +
  geom_histogram(bins = 30, color = "grey30", fill = "grey80") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),
    axis.title = element_blank()
  )


# transformation 
joined_df_tr <- joined_df %>%
  drop_na() %>%
  dplyr::select(Protein.conc.ug.ml, cell.cm, chl.cell, Fv.Fm, Sigma, Pmax.e.s, p)  %>%
  mutate(
    cell.cm = sqrt(cell.cm),
    Protein.conc.ug.ml = sqrt(Protein.conc.ug.ml),
    Pmax.e.s = log10(Pmax.e.s),
    chl.cell = log10(chl.cell)
  )

joined_df_tr %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  ) -> dat_long_tr

ggplot(dat_long_tr, aes(x = value)) +
  geom_histogram(bins = 30, color = "grey30", fill = "grey80") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),
    axis.title = element_blank()
  )


#### PCA of July and November ####
# preparing data
joined_df <- full_join(physio, fire, by = c("coral.number", "Treatment", "Coral"))

joined_df  %>%
  drop_na() %>%
  dplyr::select(Treatment, Protein.conc.ug.ml, cell.cm, chl.cell, Fv.Fm, Sigma, Pmax.e.s, p) %>%
  mutate(
    cell.cm = sqrt(cell.cm),
    Protein.conc.ug.ml = sqrt(Protein.conc.ug.ml),
    Pmax.e.s = log10(Pmax.e.s),
    chl.cell = log10(chl.cell)
  ) -> pca.physio
str(pca.physio)

# Keep only numeric variables for PCA
param_data <- pca.physio[, -which(colnames(pca.physio) %in% c("Treatment"))]
# pca
pca_result <- rda(param_data, scale = TRUE)   # RDA with no constraints = PCA

# Percent variance explained
eig_vals <- eigenvals(pca_result)
var_explained <- round(100 * eig_vals / sum(eig_vals), 1)

#  Extract sample scores 
pca_scores <- as.data.frame(scores(pca_result, display = "sites", choices = 1:3))
pca_scores$Treatment <- pca.physio$Treatment

# Extract PCA loadings (physiological variables) 
loadings_physio <- as.data.frame(scores(pca_result, display = "species", choices = 1:3))
loadings_physio$Parameter <- rownames(loadings_physio)

# Fit environmental vectors (temperature & PAR only) 
temp_par <- read.csv("../loggers/temp.par.nmds.all.csv", stringsAsFactors = F)

physio.partemp <- left_join(joined_df, temp_par, by = c("coral.number", "Treatment", "Coral"))
physio.partemp %>%
  drop_na() %>%
  dplyr::select(Treatment, Coral, temperature, PAR) -> pca.partemp

env_fit_par <- envfit(pca_result, pca.partemp[, c("temperature", "PAR")], perm = 999, choices = 1:3)
vectors_par <- as.data.frame(scores(env_fit_par, "vectors", choices = 1:3))
vectors_par$Parameter <- rownames(vectors_par)
vectors_par$R2 <- round(env_fit_par$vectors$r, 2)
vectors_par$pval <- round(env_fit_par$vectors$pvals, 3)

loadings_physio <- loadings_physio %>%
  mutate(
    Label = dplyr::recode(
      Parameter,
      Protein.conc.ug.ml = "Protein",
      cell.cm = "Cell/cm",
      chl.cell = "Chlorophyll/cell",
      Fv.Fm = "Fv/Fm",
      Sigma = "σPSII",
      Pmax.e.s = "Pmax",
      p = "p"
    )
  )

# --- Plot PCA ---
pca_scores$PC1 <- -pca_scores$PC1
loadings_physio$PC1 <- -loadings_physio$PC1
vectors_par$PC1 <- -vectors_par$PC1

pt_to_mm <- function(pt) pt / 2.845

treatment_colors <- c(
  "S" = "#b495ff",
  "SS" = "#e16cb2",
  "SD" = "#4d68a4",
  "D" = "#9aced3",
  "DD" = "#58b799",
  "DS" = "#f97642"
)

pcap <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 1) +
  #stat_ellipse(level = 0.95) +
  # Physiological loadings
  geom_segment(data = loadings_physio, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.1, "cm")), color = "black", 
               linewidth = 0.3) +
  geom_text(data = loadings_physio, aes(x = PC1, y = PC2, label = Label),
            vjust = -0.5, hjust = 0.5, size = pt_to_mm(7), color = "black") +
  # Environmental vectors (fitted)
  geom_segment(data = vectors_par, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.1, "cm")), color = "red", 
               linewidth = 0.3) +
  geom_text(data = vectors_par, aes(x = PC1, y = PC2, label = Parameter),
            vjust = -0.5, hjust = 0.5, size = pt_to_mm(7), color = "red") +
  scale_color_manual(
    values = treatment_colors,
    name = "Treatment"
  )+
  labs(
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal()

pcap
# ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/pca_nov2_revision.jpg", pcap, width = 8, height = 6)
# ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/pca_nov2_revision.jpg", pcap, width = 8, height = 6)

#### PCA of July ####
# prepare data
names(growth) <- c("growth", "coral.number", "Treatment", "Coral", "band", "rest", "origin", "sample_id", "colony_id")

physio.growth <- left_join(joined_df, growth,  by = c("coral.number", "Treatment", "Coral"))

# join everything into one table first
pca.full <- joined_df %>%
  left_join(growth, by = c("coral.number", "Treatment", "Coral")) %>%
  left_join(temp_par, by = c("coral.number", "Treatment", "Coral")) %>%
  dplyr::select(
    Treatment,
    Protein.conc.ug.ml,
    cell.cm,
    chl.cell,
    Fv.Fm,
    Sigma,
    Pmax.e.s,
    p,
    growth,
    temperature,
    PAR
  ) %>%
  drop_na()

# transformed physiological data for PCA
pca.physio <- pca.full %>%
  mutate(
    cell.cm = sqrt(cell.cm),
    Protein.conc.ug.ml = sqrt(Protein.conc.ug.ml),
    Pmax.e.s = log10(Pmax.e.s),
    chl.cell = log10(chl.cell)
  )
pca.physio <- pca.physio %>%
  mutate(Treatment = factor(Treatment, levels = c("SS", "SD", "DD", "DS")))
# PCA variables only
param_data <- pca.physio %>%
  dplyr::select(
    Protein.conc.ug.ml,
    cell.cm,
    chl.cell,
    Fv.Fm,
    Sigma,
    Pmax.e.s,
    p,
    growth
  )

# PCA
pca_result <- rda(param_data, scale = TRUE)

# percent variance explained
eig_vals <- eigenvals(pca_result)
var_explained <- round(100 * eig_vals / sum(eig_vals), 1)

# sample scores
pca_scores <- as.data.frame(scores(pca_result, display = "sites", choices = 1:3))
pca_scores$Treatment <- pca.physio$Treatment

# physiological loadings
loadings_physio <- as.data.frame(scores(pca_result, display = "species", choices = 1:3))
loadings_physio$Parameter <- rownames(loadings_physio)

loadings_physio <- loadings_physio %>%
  mutate(
    Label = dplyr::recode(
      Parameter,
      Protein.conc.ug.ml = "Protein",
      cell.cm = "Cell/cm",
      chl.cell = "Chlorophyll/cell",
      Fv.Fm = "Fv/Fm",
      Sigma = "σPSII",
      Pmax.e.s = "Pmax",
      p = "p",
      growth = "linear extension"
    )
  )

# environmental variables from the SAME rows
env_data <- pca.physio %>%
  dplyr::select(temperature, PAR)

env_fit_par <- envfit(
  pca_result,
  env_data,
  perm = 999,
  choices = 1:3
)

vectors_par <- as.data.frame(scores(env_fit_par, "vectors", choices = 1:3))
vectors_par$Parameter <- rownames(vectors_par)
vectors_par$R2 <- round(env_fit_par$vectors$r, 2)
vectors_par$pval <- round(env_fit_par$vectors$pvals, 3)
#  Plot PCA 

#flip
pca_scores$PC2 <- -pca_scores$PC2
loadings_physio$PC2 <- -loadings_physio$PC2
vectors_par$PC2 <- -vectors_par$PC2

pca_scores$PC1 <- -pca_scores$PC1
loadings_physio$PC1 <- -loadings_physio$PC1
vectors_par$PC1 <- -vectors_par$PC1

pcap_july <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 1) +
  #stat_ellipse(level = 0.95) +
  # Physiological loadings
  geom_segment(data = loadings_physio, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.1, "cm")), color = "black", 
               linewidth = 0.3) +
  geom_text(data = loadings_physio, aes(x = PC1, y = PC2, label = Label),
            vjust = -0.5, hjust = 0.5, size = pt_to_mm(7), color = "black") +
  # Environmental vectors (fitted)
  geom_segment(data = vectors_par, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.1, "cm")), color = "red", 
               linewidth = 0.3) +
  geom_text(data = vectors_par, aes(x = PC1, y = PC2, label = Parameter),
            vjust = -0.5, hjust = 0.5, size = pt_to_mm(7), color = "red") +
  scale_color_manual(
    values = treatment_colors,
    name = "Treatment"
  )+
  labs(
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal()

pcap_july
# ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/pca2_revision.jpg", pcap, width = 8, height = 6)
# ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/pca2_revision.pdf", pcap, width = 8, height = 6)
# 

#  vectors significance 
vectors_par

# pc significance
eig_vals <- eigenvals(pca_result)
var_expl <- 100 * eig_vals / sum(eig_vals)
scree_df <- data.frame(PC = paste0("PC", 1:length(var_expl)),
                       Variance = var_expl)

ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = round(Variance,1)), vjust=-0.5) +
  labs(title = "Scree plot of PCA", y = "Variance explained (%)") +
  theme_minimal()

#### loadings ####
# run for each pca result
# loadings
infl <- data.frame(
  Variable = rownames(loadings_physio),
  PC1 = loadings_physio[,1],
  PC2 = loadings_physio[,2]
)
infl[order(infl$PC1), ]

infl_plot <- infl %>%
  arrange(PC1) %>%
  mutate(Variable = factor(Variable, levels = Variable))

ggplot(infl_plot, aes(x = Variable, y = PC1)) +
  geom_col(fill = "grey30") +
  coord_flip() +
  theme_minimal() +
  labs(
    x = NULL,
    y = "Loading on PC1",
  #  title = "Physiological variable loadings on PC1"
  )


infl_long <- infl %>%
  pivot_longer(
    cols = c(PC1, PC2),
    names_to = "PC",
    values_to = "loading") %>%
  mutate(
        Label = dplyr::recode(
          Variable,
          Protein.conc.ug.ml = "Protein",
          cell.cm = "Cell/cm",
          chl.cell = "Chlorophyll/cell",
          Fv.Fm = "Fv/Fm",
          Sigma = "σPSII",
          Pmax.e.s = "Pmax",
          p = "p",
          growth = "linear extension"
        )
  )

infl %>%
  mutate(contrib = sqrt(PC1^2 + PC2^2)) %>%
  mutate(
    Label = dplyr::recode(
      Variable,
      Protein.conc.ug.ml = "Protein",
      cell.cm = "Cell/cm",
      chl.cell = "Chlorophyll/cell",
      Fv.Fm = "Fv/Fm",
      Sigma = "σPSII",
      Pmax.e.s = "Pmax",
      p = "p",
      growth = "linear extension"
    )
  )%>%
  arrange(contrib) %>%
  mutate(Label = factor(Label, levels = Label)) %>%
  pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "loading") %>%
  ggplot(aes(x = loading, y = Label, color = PC)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = loading, yend = Label), linewidth = 0.4) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("PC1" = "darkorange3", "PC2" = "darkorchid3")) +
  theme_minimal() +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 10)
  ) +
  labs(
    x = "PCA loadings",
    y = NULL,
    color = NULL,
   # title = "Correlation between traits and PC1-PC2"
  ) -> load

load

#### loadings July ####
# loadings
infl <- data.frame(
  Variable = rownames(loadings_physio),
  PC1 = loadings_physio[,1],
  PC2 = loadings_physio[,2]
)
infl[order(infl$PC1), ]

infl_plot <- infl %>%
  arrange(PC1) %>%
  mutate(Variable = factor(Variable, levels = Variable))

ggplot(infl_plot, aes(x = Variable, y = PC1)) +
  geom_col(fill = "grey30") +
  coord_flip() +
  theme_minimal() +
  labs(
    x = NULL,
    y = "Loading on PC1",
    #  title = "Physiological variable loadings on PC1"
  )


infl_long <- infl %>%
  pivot_longer(
    cols = c(PC1, PC2),
    names_to = "PC",
    values_to = "loading") %>%
  mutate(
    Label = dplyr::recode(
      Variable,
      Protein.conc.ug.ml = "Protein",
      cell.cm = "Cell/cm",
      chl.cell = "Chlorophyll/cell",
      Fv.Fm = "Fv/Fm",
      Sigma = "σPSII",
      Pmax.e.s = "Pmax",
      p = "p",
      growth = "linear extension"
    )
  )

infl %>%
  mutate(contrib = sqrt(PC1^2 + PC2^2)) %>%
  mutate(
    Label = dplyr::recode(
      Variable,
      Protein.conc.ug.ml = "Protein",
      cell.cm = "Cell/cm",
      chl.cell = "Chlorophyll/cell",
      Fv.Fm = "Fv/Fm",
      Sigma = "σPSII",
      Pmax.e.s = "Pmax",
      p = "p",
      growth = "linear extension"
    )
  )%>%
  arrange(contrib) %>%
  mutate(Label = factor(Label, levels = Label)) %>%
  pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "loading") %>%
  ggplot(aes(x = loading, y = Label, color = PC)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = loading, yend = Label), linewidth = 0.4) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("PC1" = "darkorange3", "PC2" = "darkorchid3")) +
  theme_minimal() +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 5)
  ) +
  labs(
    x = "PCA loadings",
    y = NULL,
    color = NULL,
    # title = "Correlation between traits and PC1-PC2"
  ) -> load_july

load_july

#### pca combining plots ####
fig_theme_sup <- theme(
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 7),
  text = element_text(size = 7)
)

pcap2 <- pcap + fig_theme_sup
pcap_july2 <- pcap_july + fig_theme_sup

combined <- pcap2 + pcap_july2
combined <- combined + 
  plot_annotation(tag_levels = "A")& 
  theme(plot.tag = element_text(face = "bold", size = 9)) 
combined

ggsave(
  "/home/gospozha/haifa/cayman/manuscript/revision/figs/combined_fig4.pdf",
  plot = combined,
  width = 18.4,
  height = 9,
  units = "cm",
  device = "pdf",
  useDingbats = FALSE
)
ggsave(
  "/home/gospozha/haifa/cayman/manuscript/revision/figs/combined_fig4.jpg",
  plot = combined,
  width = 18.4,
  height = 9,
  units = "cm",
  device = "jpg"
)
# # for PCAs
# ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/combined_pca.revision.jpg", 
#        combined, width = 12, height = 6)
# ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/combined_pca.revision.pdf", 
#        combined, width = 12, height = 6)
# 
# # for July
# ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/combined_pca.jul.lab.jpg", 
#        combined, width = 12, height = 6)
# ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/combined_pca.jul.lab.pdf", 
#        combined, width = 12, height = 6)


load2 <- load + fig_theme_sup
load_july2 <- load_july + fig_theme_sup
combined <- load2 + load_july2
combined <- combined + 
  plot_annotation(tag_levels = "A")& 
  theme(plot.tag = element_text(face = "bold")) 
combined

ggsave(
  "/home/gospozha/haifa/cayman/manuscript/revision/figs/combined_figS7.pdf",
  plot = combined,
  width = 18.4,
  height = 9,
  units = "cm",
  device = "pdf",
  useDingbats = FALSE
)
ggsave(
  "/home/gospozha/haifa/cayman/manuscript/revision/figs/combined_figS7.jpg",
  plot = combined,
  width = 18.4,
  height = 9,
  units = "cm",
  device = "jpg"
)
# # for PCAs
# ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/combined_load.revision.jpg", 
#        combined, width = 12, height = 6)
# ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/combined_load.revision.pdf", 
#        combined, width = 12, height = 6)
# 

#### permanova Nov and Jul ####
# run for each pca result
# nov and july
partemp <- left_join(joined_df, temp_par, by = c("coral.number", "Treatment", "Coral"))
#july
#partemp <- left_join(physio.growth, temp_par, by = c("coral.number", "Treatment", "Coral"))
partemp %>%
  drop_na() %>%
  dplyr::select(Treatment, temperature, PAR, Protein.conc.ug.ml, cell.cm, chl.cell, Fv.Fm, Sigma, Pmax.e.s, p) %>%
  #dplyr::select(Treatment, temperature, PAR, growth, Protein.conc.ug.ml, cell.cm, chl.cell, Fv.Fm, Sigma, Pmax.e.s, p) %>%
  mutate(
    cell.cm = sqrt(cell.cm),
    Protein.conc.ug.ml = sqrt(Protein.conc.ug.ml),
    Pmax.e.s = log10(Pmax.e.s),
    chl.cell = log10(chl.cell)
  )  -> perm.partemp
str(perm.partemp)

perm.sc <- perm.partemp[ , -c(1,2,3)] %>% scale()

permanova_treatment <- adonis2(
  perm.sc ~ Treatment,
  data = perm.partemp,
  method = "euc",
  perm=999
)
permanova_treatment

pairwise_res <- pairwise.adonis(
  x =perm.sc,
  factors = perm.partemp$Treatment,
  sim.method="euclidean",
  perm = 999
)

pairwise_res



# 
# # numeric matrix exactly as used for PCA
# physio_mat <- pca.physio %>%
#   dplyr::select(-Treatment) %>%
#   scale()
# 
# # Euclidean distance (correct for scaled continuous traits)
# dist_physio <- dist(physio_mat, method = "euclidean")
# 
# permanova_treatment <- adonis2(
#   dist_physio ~ Treatment,
#   data = pca.physio,
#   permutations = 999
# )
# 
# permanova_treatment
# 
# 
# pairwise_res <- pairwise.adonis(
#   x = dist_physio,
#   factors = pca.physio$Treatment,
#   perm = 999
# )
# 
# pairwise_res

dist_mat <- dist(perm.sc, method = "euclidean")

disp <- betadisper(dist_mat, perm.partemp$Treatment)

anova(disp)
permutest(disp, permutations = 999)


#### permanova for July ####
partemp <- left_join(
  physio.growth,
  temp_par,
  by = c("coral.number", "Treatment", "Coral")
)

perm.partemp <- partemp %>%
  drop_na() %>%
  dplyr::select(
    Treatment, temperature, PAR,
    growth, Protein.conc.ug.ml, cell.cm, chl.cell,
    Fv.Fm, Sigma, Pmax.e.s, p
  ) %>%
  filter(Treatment %in% c("SS", "SD", "DD", "DS")) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("SS", "SD", "DD", "DS")),
    
    origin = case_when(
      Treatment %in% c("SS", "SD") ~ "Shallow",
      Treatment %in% c("DD", "DS") ~ "Deep"
    ),
    
    transplant_depth = case_when(
      Treatment %in% c("SS", "DS") ~ "Shallow",
      Treatment %in% c("SD", "DD") ~ "Deep"
    ),
    
    origin = factor(origin, levels = c("Shallow", "Deep")),
    transplant_depth = factor(transplant_depth, levels = c("Shallow", "Deep")),
    
    cell.cm = sqrt(cell.cm),
    Protein.conc.ug.ml = sqrt(Protein.conc.ug.ml),
    Pmax.e.s = log10(Pmax.e.s),
    chl.cell = log10(chl.cell)
  )

# response matrix: physiological traits only
perm.sc <- perm.partemp %>%
  dplyr::select(
    growth, Protein.conc.ug.ml, cell.cm, chl.cell,
    Fv.Fm, Sigma, Pmax.e.s, p
  ) %>%
  scale()

# PERMANOVA: origin × transplant depth
permanova_origin_transplant <- adonis2(
  perm.sc ~ origin * transplant_depth,
  data = perm.partemp,
  method = "euclidean",
  by = "term",
  permutations = 999
)

permanova_origin_transplant


pairwise_res_treatment <- pairwise.adonis(
  x = perm.sc,
  factors = perm.partemp$Treatment,
  sim.method = "euclidean",
  perm = 999
)

pairwise_res_treatment


dist_sc <- dist(perm.sc, method = "euclidean")

disp_treatment <- betadisper(dist_sc, perm.partemp$Treatment)
permutest(disp_treatment, permutations = 999)

perm.partemp$group <- interaction(
  perm.partemp$origin,
  perm.partemp$transplant_depth,
  sep = "_"
)

disp_group <- betadisper(dist_sc, perm.partemp$group)
permutest(disp_group, permutations = 999)
#### survival  ####

surv <- read.csv("../NSF_transplant_survival.csv", stringsAsFactors = FALSE)

surv <- surv %>%
  mutate(
    Treatment = factor(Treatment, levels = c("SS", "SD", "DD", "DS")),
    Site      = factor(Site),
    Status    = factor(Status),
    Dead      = ifelse(Status %in% c("D"), 1L, 0L)
  )

print(table(surv$Site, surv$Treatment, surv$Dead))

## Quick global tests (within site) + effect sizes
tab_site <- split(surv, surv$Site)

chisq_or_fisher <- function(df) {
  tt <- table(df$Treatment, df$Dead)
  # if any expected cell < 5, use Fisher (safe default)
  if (any(chisq.test(tt)$expected < 5)) {
    out <- fisher.test(tt)
    list(test="Fisher", p.value=out$p.value, odds.ratio=unname(out$estimate))
  } else {
    out <- chisq.test(tt)
    list(test="Chi-square", p.value=out$p.value)
  }
}

by_site_tests <- lapply(tab_site, chisq_or_fisher)
print(by_site_tests)
## Quick global tests (within site) + effect sizes
tab_site <- split(surv, surv$Site)

chisq_or_fisher <- function(df) {
  tt <- table(df$Treatment, df$Dead)
  # if any expected cell < 5, use Fisher (safe default)
  if (any(chisq.test(tt)$expected < 5)) {
    out <- fisher.test(tt)
    list(test="Fisher", p.value=out$p.value, odds.ratio=unname(out$estimate))
  } else {
    out <- chisq.test(tt)
    list(test="Chi-square", p.value=out$p.value)
  }
}

by_site_tests <- lapply(tab_site, chisq_or_fisher)
print(by_site_tests)
#$`Coral City`$p.value
#[1] 0.1667579

#$`Martha's Finyard`$p.value
#[1] 0.03311941

## Logistic regression with Treatment * Site (binomial)
m_glm <- glm(Dead ~ Treatment * Site, data = surv, family = binomial())
summary(m_glm)

## Type-II/III tests for main effects + interaction
car::Anova(m_glm, type = 2)  # robust default; switch to type=3 if you prefer

m_glm2 <- glm(Dead ~ Treatment + Site, data = surv, family = binomial())
summary(m_glm2)

## Type-II/III tests for main effects + interaction
car::Anova(m_glm2, type = 2)  # robust default; switch to type=3 if you prefer

AIC(m_glm, m_glm2)

# exclude interaction
emm_trt <- emmeans(m_glm2, ~ Treatment)

pairs_trt <- pairs(emm_trt, adjust = "tukey")
pairs_trt

emm_site <- emmeans(m_glm2, ~ Site)

pairs_site <- pairs(emm_site, adjust = "tukey")
pairs_site


# Prepare data with percentages
plot_data <- surv %>%
  dplyr::count(Site, Treatment, Dead, name = "N") %>%
  dplyr::group_by(Site, Treatment) %>%
  dplyr::mutate(Percent = 100 * N / sum(N),
                Status = ifelse(Dead == 1, "Dead", "Alive"))

# Patterns 
fill_patterns <- c(
  "S"       = "none",
  "SS"    = "circle",
  "SD"    = "stripe",   # pattern like 10→10
  "D"       = "none",
  "DD"    = "stripe",
  "DS"    = "circle" # pattern like 40→40
)

# Basic stacked barplot
plot_data <- plot_data %>%
  mutate(
    fill_color = case_when(
      Status == "Dead" ~ "gray30",
      Treatment == "SS" & Status == "Alive" ~ "#f75f55",
      Treatment == "DS" & Status == "Alive" ~ "#00A9FF",
      Treatment == "DD" & Status == "Alive" ~ "#00A9FF",
      Treatment == "SD" & Status == "Alive" ~ "#f75f55",
    )
  )

s <- ggplot(plot_data, aes(x = Treatment, y = Percent, 
                           fill = fill_color, pattern = Treatment)) +
  geom_bar_pattern(stat = "identity",
                   aes(group = Status),
                   pattern_density = 0.08,
                   pattern_spacing = 0.02,
                   pattern_fill = "grey",
                   color = "gray30") +
  facet_wrap(~ Site, ncol = 2) +
  scale_fill_identity() +  # Use raw hex colors, no scale mapping
  scale_pattern_manual(values = fill_patterns) +
  ylab("") +
  xlab("") +
  #ggtitle("Survival rate (Alive/Dead)") +
  #scale_x_discrete(labels = c("10→10", "10→40", "40→40", "40→10")) +
  theme_minimal() +
  geom_text(
    aes(label = paste0(round(Percent, 1), "%"), group = Status),
    position = position_stack(vjust = 0.5),
    size = 3.7, color = "white"
  ) +
  guides(pattern = "none")  
s
ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/survival_horiz2.jpg", s, width = 8, height = 5)
ggsave("~/haifa/cayman/P.astreoides_physiology/github3/5/survival_horiz2.pdf", s, width = 8, height = 5)




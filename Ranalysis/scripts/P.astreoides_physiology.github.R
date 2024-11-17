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

# working directory
setwd("~/haifa/cayman/P.astreoides_physiology/")

# opening file, checking, replacing characters to vectors when needed
physio <- read.csv('Physiology_P.astreoides.csv', stringsAsFactors = F)
str(physio)
View(physio)
physio$Treatment = factor(physio$Treatment, levels = c("S.T0", "D.T0", "SS", "SD", "DD", "DS"))
physio$Coral = factor(physio$Coral, levels = c("MF", "CC"))
physio$Month = factor(physio$Month, levels = c("Nov", "July"))

# adding the same theme to each plot
mytheme = theme_bw()+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 8), axis.title = element_text(size = 9))

# analysis

### Protein.cm ###

# making basic ggplot object 
p.protein = ggplot(physio, aes(y = Protein.conc.ug.ml, x = Treatment)) +
  labs(y= "Protein concentration (ug/ml)")+
  mytheme +
  guides(fill = FALSE)

# drawing boxplots of each treatment for each location
p.protein + 
  geom_boxplot(outlier.shape = NA, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.4)  +
  geom_jitter(position = position_jitter(width = .25), size = 0.5) +
  facet_wrap(physio$Coral)

# checking for normality and homoscedasticity
shapiro.test(physio$Protein.conc.ug.ml)  # data is not normal (p.value < 0.5)
leveneTest(Protein.conc.ug.ml~Treatment, d=physio) # heteroscedasticity of variance (p.value < 0.5)

# trying different data transformations
shapiro.test(log(physio$Protein.conc.ug.ml))
shapiro.test(sqrt(physio$Protein.conc.ug.ml)) # data is normal after sqrt transformation
leveneTest(sqrt(Protein.conc.ug.ml)~Treatment, d=physio) # variance is almost ok, we need to check QQ plots

# performing anova + Tukey post-hoc 

# starting with the full model (includes treatment, location and their interaction)
model <- lm(data = physio, sqrt(Protein.conc.ug.ml) ~ Treatment * Coral)
# backward stepwise selection: a comparison of models based on the AIC criterion
drop1(model) # model without interaction is better
# dropping the interaction
model1 <- update(model, .~. - Treatment:Coral)  
drop1(model1) # model without Coral is better
# dropping the location
model2 <- update(model1, .~. - Coral)
# checking summary and plots
summary(model2)
plot(model2) # Q-Q plots - checking for heteroscedasticity of residuals - look very good
# data is ok for anova
anova(model2)

# posthoc test - pairwise comparisons (Tukey test)
posthoc <- TukeyHSD(aov(model2))
posthoc$Treatment # the table of pairwise comparisons
# assigning letters to groups ("compact letter display")
letters <- multcompLetters4(aov(model2), posthoc) 
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
p.protein + 
  geom_boxplot(outlier.size = 0.5, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.4)  + 
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -1.25, vjust = -0.8, fontface = "italic")


### cell/cm ###

# plots
p.cellcm = ggplot(physio, aes(y = cell.cm, x = Treatment)) +
  labs(y = ~cells ~x10^5 ~cm^-2)+
  mytheme +
  guides(fill = FALSE)

p.cellcm  +
  geom_boxplot(outlier.shape = NA, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.5) +
  facet_wrap(physio$Coral)

# normality, homoscedasticity
shapiro.test(physio$cell.cm) # not normal  
leveneTest(cell.cm~Treatment,d=physio) # variance is ok

#transformation
shapiro.test(sqrt(physio$cell.cm)) # normal 
leveneTest(sqrt(cell.cm)~Treatment,d=physio) # variance is ok

# anova + Tukey post-hoc 

# starting with full model
model <- lm(data = physio, sqrt(cell.cm) ~ Treatment * Coral)
# backward stepwise selection
drop1(model) 
# dropping the interaction because it's not significant
model1 <- update(model, .~. - Treatment:Coral)
drop1(model1) # Coral is also significant, we can't drop it. But the difference in AIC is very small.
# It means that the location explains a little bit of variance. 

summary(model1)  -
plot(model1) # Q-Q plot is ok
anova(model1) # Coral is not significant factor 

# Tukey posthoc test
posthoc <- TukeyHSD(aov(model1)) 
posthoc$Coral # shows that CC is the same as MF
posthoc$Treatment # shows pairwise comparisons taking into account Coral in the model
letters <- multcompLetters4(aov(model1), posthoc)
letters.df <- data.frame(letters$Treatment$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Treatment <- rownames(letters.df) 
placement <- physio %>%
  group_by(Treatment) %>%
  summarise(quantile(cell.cm, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement)

# boxplot with letters
p.cellcm + 
  geom_boxplot(outlier.size = 0.5, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.4)  +
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -1.25, vjust = -0.8, fontface = "italic")


###chl/cell###

# plots
p.chlcell = ggplot(physio, aes(y = chl.cell, x = Treatment)) +
  labs(y= ~chlorophyll[a] ~x10^-7 ~cell^-1)+
  mytheme +
  guides(fill = FALSE)+
  # making y axis shorter because of outliers
  scale_y_continuous(limits = c(0,50))

p.chlcell   +
  geom_boxplot(outlier.shape = NA, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.2)  +
  geom_jitter(position = position_jitter(width = .25), size = 0.5) +
  facet_wrap(physio$Coral)

# normality, homoscedasticity
shapiro.test(physio$chl.cell) # data is not normal (p.value < 0.5)
leveneTest(chl.cell~Treatment,d=physio)  # heteroscedasticity of variance (p.value < 0.5)

# transforming the data (log, log10, sqrt, ^(1/3), exp,....)
shapiro.test(sqrt(physio$chl.cell)) # data is not normal (p.value < 0.5)
leveneTest(sqrt(chl.cell)~Treatment,d=physio)  # heteroscedasticity of variance (p.value < 0.5)
# transformations are not working

# Kruskal-Wallis test for nonparametric data

# cannot add in the model the second factor (Coral) simultaneously
# checking location first 
kruskal.test(chl.cell ~ Coral, data = physio)
# location is significant

# performing separate tests for each location

# MF
kruskal.test(chl.cell ~ Treatment, data = physio[physio$Coral=='MF',]) #  significant

# posthoc: Dunn's Test with Bonferroni correction for p-values
dunn.res <- dunnTest(chl.cell ~ Treatment,
                     data=physio[physio$Coral=='MF',],
                     method="bonferroni")
diff <- dunn.res$res$P.adj < 0.05
Names <- gsub(" ", "", dunn.res$res$Comparison)
names(diff) <- Names
# compact letter display
letters <- multcompLetters(diff)
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter" 
letters.df$Treatment <- rownames(letters.df) 
placement <- physio[physio$Coral=='MF',] %>%
  group_by(Treatment) %>%
  summarise(quantile(chl.cell, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 

# boxplot with letters

ggplot(physio[physio$Coral=='MF',], aes(y = chl.cell, x = Treatment)) +
  geom_boxplot(outlier.size = 0.5, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  labs(y= ~chlorophyll[a] ~x10^-7 ~cell^-1)+
  mytheme +
  guides(fill = FALSE)+
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -1.25, vjust = -0.8, fontface = "italic")


# CC
kruskal.test(chl.cell ~ Treatment, data = physio[physio$Coral=='CC',]) #  significant

# posthoc: Dunn's Test with Bonferroni correction for p-values
dunn.res <- dunnTest(chl.cell ~ Treatment,
                     data=physio[physio$Coral=='CC',],
                     method="bonferroni")
diff <- dunn.res$res$P.adj < 0.05
Names <- gsub(" ", "", dunn.res$res$Comparison)
names(diff) <- Names
# compact letter display
letters <- multcompLetters(diff)
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Treatment <- rownames(letters.df) 
placement <- physio[physio$Coral=='CC',] %>% 
  group_by(Treatment) %>%
  summarise(quantile(chl.cell, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 
letters.df[3, 3] <- 13 # changing this position because of the outlier

# boxplot with letters

ggplot(physio[physio$Coral=='CC',], aes(y = chl.cell, x = Treatment)) +
  geom_boxplot(outlier.size = 0.5, aes(fill = Treatment), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  labs(y= ~chlorophyll[a] ~x10^-7 ~cell^-1)+
  mytheme +
  guides(fill = FALSE)+
  # making y axis shorter because of outliers
  scale_y_continuous(limits = c(0,50))+ 
  geom_text(data = letters.df, aes(x = Treatment, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")


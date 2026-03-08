# Adapted from Carpenter et al 2022
# Light and photoacclimatization drive distinct differences between shallow and mesophotic coral communities



fire.data <- read.csv("/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Cayman-translocation/Ranalysis/Data/FIRe_tranplantation_exsitu_Nov2022_July2023.csv")

#FIRe Data

library(tidyverse)
install.packages("lubridate")
library(broom)
library(lubridate)
library(cowplot)
library(tidyverse)
library(ggtext)

fire.data$depth <- as_factor(fire.data$depth)


is.na(fire.data$depth) #Line 421 - line 425
fire.data <- na.omit(fire.data)

fvfm <- ggplot(data = fire.data,aes(x= depth, y=Fv/Fm, fill=Species))+
  geom_boxplot() +
  theme_set(theme_cowplot()) +
  labs(title="Quantum yield of photochemistry in PSII", y="Fv’/Fm’") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  #scale_fill_manual(values=c("#4393C3","#F4A582")) +
  scale_x_discrete(limits = rev) +
  coord_flip()

fvfm

sigma <- ggplot(data = fire.data,aes(x= depth, y=Sigma, fill=Species))+
  geom_boxplot() +
  theme_set(theme_cowplot()) +
  labs(title="Functional absorption cross-section of PSII", y= "σPSII’(A2)") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  #  scale_fill_manual(values=c("#4393C3","#F4A582")) +
  scale_x_discrete(limits = rev) +
  coord_flip()
sigma

pmax <- ggplot(data = fire.data,aes(x= depth, y=Pmax.e.s, fill=Species))+
  geom_boxplot() +
  theme_set(theme_cowplot()) +
  labs(title="Maximum photosynthetic rate ", y="Pmax (electron s-1 PSII-1)") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  #  scale_fill_manual(values=c("#4393C3","#F4A582")) +
  scale_x_discrete(limits = rev) +
  coord_flip()
pmax

p <- ggplot(data = fire.data,aes(x= depth, y=p, fill=Species))+
  geom_boxplot() +
  theme_set(theme_cowplot()) +
  labs(title="connectivity parameter ", y="p") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  #scale_fill_manual(values=c("#4393C3","#F4A582")) +
  scale_x_discrete(limits = rev) +
  coord_flip()
p

p_legend <- ggplot(data = fire.data,aes(x= depth, y=p, fill=Species))+
  geom_boxplot() +
  theme_set(theme_cowplot()) +
  labs(title="p") +
  #  scale_fill_manual(values=c("#4393C3","#F4A582")) +
  coord_flip()
p_legend

nolegend <- plot_grid(fvfm,sigma,pmax, p, labels = c('a','b','c','d'))

legend <- get_legend(p_legend)
fire_plot <- plot_grid(nolegend,legend, rel_widths = c(1,0.25))
fire_plot



#FIRe STATS
```{r}
library("FSA")
library("rcompanion")

#Kruskal-Wallis Test- data is not normally distributed
fire.data$Depth_group <- as_factor(fire.data$depth)
levels(fire.data$Depth_group)
fire.past <- subset(fire.data, Species == "Porites astreoides")


#FvFm
#First, look at differences between species
KW_fvfm_s <- kruskal.test(Fv/Fm ~ Species, data = fire.data) #null = there is no difference in value
KW_fvfm_s
#Then, look at differences between depts for each species
KW_fvfm_past <- kruskal.test(Fv/Fm ~ depth, data = fire.past)
KW_fvfm_past
DT_fvfm_past <- dunnTest(fire.past$"Fv/Fm" ~ fire.past$depth, method = "bonferroni")
DT_fvfm_past
fvfm_past_let <- DT_fvfm_past$res



#Then, look at differences between depts for each species
KW_sigma_past <- kruskal.test(Sigma ~ depth, data = fire.past)
KW_sigma_past #no sig difference, so no post hoc




# sigma past

KW_sigma_past <- kruskal.test(Sigma ~ depth, data = fire.past)
KW_sigma_past
DT_sigma_mcav <- dunnTest(fire.mcav$Sigma ~ fire.mcav$depth, method = "bonferroni")
DT_sigma_mcav

#Then, look at differences between depts for each species
KW_p_past <- kruskal.test(p ~ depth, data = fire.past)
KW_p_past
DT_p_past <- dunnTest(fire.past$p ~ fire.past$depth, method = "bonferroni")
DT_p_past

KW_p_mcav <- kruskal.test(p ~ depth, data = fire.mcav)
KW_p_mcav 
DT_p_past <- dunnTest(fire.past$p ~ fire.past$depth, method = "bonferroni")
DT_p_past



#Pmax
#Then, look at differences between depts for each species
KW_pmax_past <- kruskal.test(Pmax.e.s ~ depth, data = fire.past)
KW_pmax_past
DT_pmax_past <- dunnTest(fire.past$Pmax.e.s ~ fire.past$depth, method = "bonferroni")
DT_pmax_past

KW_pmax_mcav <- kruskal.test(Pmax.e.s ~ depth, data = fire.mcav)
KW_pmax_mcav
DT_pmax_mcav <- dunnTest(fire.mcav$Pmax.e.s ~ fire.mcav$depth, method = "bonferroni")
DT_pmax_mcav
```
#LI-CORE
```{r}

library(ggplot2)

licore_data <- read.csv('LICOR_06242021.csv')
licore_data$Depth <- as.factor(licore_data$Depth)
licore <-ggplot(licore_data, aes(Depth, INPUT1, fill=Reading)) +
  geom_boxplot()
licore

licore_percent_all <- licore_data %>%
  mutate(percent_of_max = INPUT1/311.12200) #Solar noon surface = 311.12200

licore_percent_plot <- ggplot(licore_percent_all, aes(x=Depth,y=percent_of_max)) +
  geom_boxplot()
licore_percent_plot

licore_means <- licore_data %>%            #  THIS IS THE DATA SET USED TO MAKE PAR DATA IN MS
  group_by(Depth, Reading) %>%
  summarise(PAR_means = mean(INPUT1),
            std_dev = sd(INPUT1)) %>%
  ungroup()

licore_percent <- read.csv('licor_percents.csv')

summary(licore_means)

licore_aov <- aov (INPUT1 ~ Reading * Depth, data=licore_data)
anova(licore_aov)
```

#MCAV GLM 
```{r}
library(dplyr)
library(ggbiplot)
library(scales)
library(ggfortify)
library(car)

fire_pca <- read.csv('FIRePCA.csv')
bloc_len <- 5
fire_pca$coral_number <- rep(seq(1,1+nrow(fire_pca)%/%bloc_len), each=bloc_len, length.out=nrow(fire_pca))

fire_pca2 <- select(fire_pca,4,17:19,45,47:50)
pca_mcav_data1 <- subset(fire_pca2, Species == "M. cavernosa")
pca_mcav_data1 <- select(pca_mcav_data1, -1)
pca_mcav_data1 <- pca_mcav_data1%>%
  group_by(coral_number) %>%
  summarise_all(~mean(.))
pca_mcav_data2 <- select(pca_mcav_data1, 2:4, 6:8)
pca_mcav_fact <- select(pca_mcav_data1, 2:4,6:7)



pca_mcav1 <- prcomp(pca_mcav_fact, center = TRUE, scale. = TRUE)
print(pca_mcav1)

mcav_factanal1<- factanal(pca_mcav_fact, factors = 2, scores = 'regression')
mcav_factanal1
pca_mcav_plot1 <- autoplot(mcav_factanal1, data = pca_mcav_data2, colour='coral_abundance', size='coral_abundance', loadings=TRUE, loadings.label = TRUE, loadings.label.vjust = -1,loadings.label.hjust=-0.1, loadings.colour = "grey50", loadings.label.colour = "black") +
  guides(color=guide_legend(), size=guide_legend()) +
  labs(colour="Coral abundance", size = "Coral abundance",) +
  theme_classic()
pca_mcav_plot1


scale(pca_mcav_data1, center = TRUE, scale = TRUE)
fit_mcav <- glm(coral_abundance~ Fv.Fm + Sigma + Pmax.e.s. + par_percent + p, data = pca_mcav_data2)
summary(fit_mcav)
vif(fit_mcav)

hist(resid(fit_mcav))
qqnorm(resid(fit_mcav)) 
qqline(resid(fit_mcav))

```

#PAST GLM 
```{r}

fire_pca2 <- select(fire_pca,4,17:19,45,47:50)
pca_past_data1 <- subset(fire_pca2, Species == "P. astreoides")
pca_past_data1 <- select(pca_past_data1, -1)
pca_past_data1 <- pca_past_data1%>%
  group_by(coral_number) %>%
  summarise_all(~mean(.))
pca_past_data2 <- select(pca_past_data1, 2:4, 6:8)
pca_past_fact <- select(pca_past_data1, 2:4,6:7)



pca_past1 <- prcomp(pca_past_fact, center = TRUE, scale. = TRUE)
print(pca_past1)

past_factanal1<- factanal(pca_past_fact, factors = 2, scores = 'regression')
pca_past_plot1 <- autoplot(past_factanal1, data = pca_past_data2, colour='coral_abundance', size='coral_abundance', loadings=TRUE, loadings.label = TRUE, loadings.label.vjust = -1,loadings.label.hjust=-0.1, loadings.colour = "grey50", loadings.label.colour = "black") +
  scale_color_fermenter(palette = "Oranges") +
  guides(color=guide_legend(), size=guide_legend()) +
  labs(colour="Coral abundance", size = "Coral abundance",) +
  theme_classic()
pca_past_plot1



scale(pca_past_data1, center = TRUE, scale = TRUE)
fit_past <- glm(coral_abundance~ Fv.Fm + Sigma + Pmax.e.s. + par_percent + p, data = pca_past_data2)
summary(fit_past)
vif(fit_past)

hist(resid(fit_past))
qqnorm(resid(fit_past)) 
qqline(resid(fit_past))
```
#Model Plots
```{r}
library(cowplot)
pca_plots <- plot_grid(pca_mcav_plot1, pca_past_plot1,labels = c('A', 'B'))
pca_plots
ggsave("Figure5.pdf", plot=pca_plots,height = 6, width = 12, units = "in")
```

#THE MAP
```{r}
library(ggplot2)
library (sf)
library(tidyverse)
library(maps)
library(ggspatial)

library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

carib <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-100.10, -74.12), ylim = c(12.10, 30.50), expand = FALSE) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y   = element_blank(),
        axis.text.x   = element_blank())

carib

ggsave("small_map.pdf", plot=carib, height = 2, width = 2, units = "in")

# function to turn a data.frame (with the columns
# "Longitude" and "Latitude") into a spatial object
tibble_to_sf <- function(tib, crs = 4326){
  tib %>%
    st_as_sf(., coords = c("long","lat")) %>%
    st_set_crs(., crs)
}


cayman <- read_sf('shapefile/gadm36_CYM_0.shp')


sites <- read.csv('LC_coordinates.csv')
sites_sf <- sites %>%
  select(site_name:long) %>%
  tibble_to_sf()

lc_map <- ggplot() +
  geom_sf(data=cayman) +
  geom_sf(data= sites_sf, aes(fill=site_name)) +
  lims(x= c(-80.15, -79.93), y = c(19.63, 19.78)) +
  annotation_scale(location = "br", width_hint = 0.5) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.3, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Longitude (W)") +
  ylab("Latitude (N)")
lc_map

ggsave("Figure_1.pdf", plot=lc_map, height = 5, width = 6, units = "in")

```

#CORRELOGRAM
```{r}
library(tidyverse)
library(dplyr)
library(vegan)


wide_data2 <- lc_sf.data %>%
  group_by(depth, species_aggra, site, file_name) %>%
  summarise(total_species = sum (coral_area_m)) %>%
  ungroup()

wide_data3 <- subset(wide_data2, site == 'dynamite_drop')

corr.data<- wide_data3 %>%
  pivot_wider(names_from = species_aggra, values_from = total_species)

corr.data_dd = corr.data[-2]

corr.data_dd <- corr.data_dd %>% 
  mutate(transect_id = str_c(depth,file_name, sep = '_')) 

corr.data_dd[is.na(corr.data_dd)] <- 0 

corr_matrix <- corr.data_dd %>% 
  select(AAGA:SLAC) %>% #need to update if species change
  as.matrix() 

rownames(corr_matrix) <- corr.data_dd$transect_id

diversity(corr_matrix,"shannon", MARGIN = 1, base = exp(1))

#PREVIOUS MATRIX AGGREGATED ALL QUADRATS, THIS ONE KEEPS THEM AS SAMPLES
#also changing from total area to count of colonies

wide_data2_sample <- lc_sf.data %>%
  group_by(depth, species_aggra, site, file_name) %>%
  summarise(n = n())

wide_data3_sample <- subset(wide_data2_sample, site == 'dynamite_drop')
#now we have the total number of species in each picture
#need to find abundance to use in matrix for shannon diversity
#so will divide species in quadrat/ total number of species at site

write.csv(wide_data3_sample,"richness.csv")

n_distinct(wide_data3_sample$species_aggra) #34 distinct species

wide_data3_sample$abundance <- wide_data3_sample$n/34

wide_data4_sample <- wide_data3_sample %>%
  group_by(depth, file_name) %>%
  summarise(n_species = n_distinct(species_aggra))



#FROM HERE ON FOR SHANNON DIVERSITY INDEX
corr.data_sample<- wide_data3_sample %>%
  pivot_wider(names_from = species_aggra, values_from = abundance)

corr.data_sample <- corr.data_sample %>% 
  mutate(transect_id = str_c(depth,file_name, sep = '_')) 

corr.data_sample = corr.data_sample[-1:-5]

corr.data_sample[is.na(corr.data_sample)] <- 0 

write.csv(corr.data_sample, "shann_corr.csv") #took it to excel to aggregate by transect_id (r code wouldn't work and I lost patience)
#aggregate(.~ transect_id, data = corr.data_sample, sum)

shann_corr <- read.csv('shann_corr.csv')

shann_corr2 <- shann_corr[-1]
shann_diversity <- diversity(shann_corr2, index = "shannon")
shann_diversity

write.csv(shann_diversity, "shann_diversity.csv")

#corr_matrix_sample <- corr.data_sample %>% 
# select(AAGA:SLAC) #%>% #need to update if species change
#as.matrix() 

#rownames(corr_matrix_sample) <- corr.data_sample$transect_id

#shann_diversity <- diversity(corr.data_sample,index = "shannon", )
#shann_diversity

#Now, rather than shannon index I am going to look at species accumulation curves

#get richness estimators (for each sample, cumulative)
#pool.species <- poolaccum(corr_matrix_sample)
#plot all: obs richness and  estimators
#plot(pool.species)

#build the species accumulation curve & rarefaction curve (expected)
species.specaccum <- specaccum(corr_matrix_sample,method = "rarefaction")
#plot the curve with some predefined settings
plot(species.specaccum,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

#build a expected curve (randomization for boxplot comparison)
species.specaccum.rand <- specaccum(corr_matrix_sample, "random")
#plot both curves ("observed" vs "randomized")
plot(species.specaccum.rand,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(species.specaccum.rand, col="yellow", add=TRUE, pch="+")

richness = coral2 %>%
  group_by(depth, site) %>%
  summarise(richness = n_distinct(species)) 

write.csv(richness, "richness.csv")

#copied this as diversity into corr matrix in excel
#import this for correlation plot
```
#iNEXT
```{r}
library(iNEXT)
species_rarecurve <- rarecurve()
```


#Using BiodiversityR and by groups
```{r}
library(BiodiversityR)
#take wide_data_3_sample that counts the number of species and adding a variable for variable

wide_data2_area <- lc_sf.data %>%
  group_by(depth, species_aggra, site, file_name, quadrat_area_m) %>%
  summarise(n = n())

wide_data3_area <- subset(wide_data2_area, site == 'dynamite_drop')

accum_data <- wide_data3_area[]

accum <- accumcomp(wide_data3_area, y=n, factor='depth', scale='quadrat_area_m',
                   method='exact', conditioned=FALSE, plotit=FALSE)
accum

##NEW CODE 1/28 AM

library(BiodiversityR) # also loads vegan
library(ggplot2)
library(ggsci)

accum_species<- wide_data3_sample %>%
  pivot_wider(names_from = species_aggra, values_from = n)

accum_species_1 <- accum_species[-1:-2]

accum_species_1[is.na(accum_species_1)] <- 0 

accum_species_mat <- accum_species_1 %>% 
  select(AAGA:SLAC) %>% #need to update if species change
  as.matrix() 

#accum_env <-wide_data3_area[!duplicated(wide_data3_area$file_name), ]

accum_env <- wide_data3_area[c(4,1,5)]

accum.1 <- accumcomp(accum_species_mat, y=accum_env, factor = 'depth',method='exact', conditioned=FALSE, plotit=FALSE)
accum.1

accum.long1 <- accumcomp.long(accum.1, ci=NA, label.freq=4)
head(accum.long1)

plotgg1 <- ggplot(data=accum.long1, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=1) +
  #geom_point(data=subset(accum.long1, labelit==TRUE), 
  #aes(colour=Grouping), size=2) +
  #geom_ribbon(aes(colour=Grouping), alpha=0.2, show.legend=FALSE) + 
  scale_colour_npg() +
  labs(x = "Quadrats", y = "Coral Species", colour = "Depth")
plotgg1

accum.2 <- accumcomp(accum_species_mat, y=accum_env, factor='depth', scale='quadrat_area_m',method='exact', conditioned=FALSE, plotit=FALSE)

accum.2

accum.long2 <- accumcomp.long(accum.2, ci=NA, label.freq=1)
head(accum.long2)

```

#Correlation Matrix
```{r}
library(ggcorrplot)

corr_matrix2 <- read.csv('correlation_matrix.csv')
corr_matrix2 <- subset(corr_matrix2, select = -c(1,3))


corr <- round(cor(corr_matrix2), 1)
corr

p.mat <- cor_pmat(corr_matrix2)
p.mat

light_corr <- ggcorrplot(corr, hc.order = TRUE, 
                         type = "lower", 
                         lab = TRUE, 
                         lab_size = 3, outline.col = "white",
                         p.mat = p.mat,
                         title="Correlogram of Light and Benthos", 
                         colors = c("darkblue", "white", "red"),
                         ggtheme=theme_bw)
light_corr
ggsave('correlogram.pdf', plot = light_corr, height = 6, width = 6, units = "in")

```

#SMOOTH PAR FIG
```{r}
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
library(RColorBrewer)

mean_par <- licore_data %>%
  group_by(Reading, Depth) %>%
  summarise(mean = mean(INPUT1), sd = sd(INPUT1)) %>%
  ungroup()

colourCount = length(unique(mean_par$Reading))
getPalette = colorRampPalette(brewer.pal(3, "Paired"))

par_line <- ggplot(mean_par, aes(x=mean,y=Depth, group=Reading,color=Reading)) +
  geom_line() +
  geom_ribbon(aes(xmin=mean-sd,xmax=mean+sd, fill = Reading), alpha=0.3) +
  scale_colour_manual(values = c("#1f78b4","#e31a1c","#fdbf6f")) +
  scale_y_discrete(limits=rev) +
  theme_classic()
par_line
ggsave('par_line.pdf', plot = par_line, height = 4, width = 6, units = "in")

par_line2 <- ggplot(mean_par, aes(x=mean,y=Depth, group=Reading)) +
  geom_line(aes(colour = factor(Reading))) +
  scale_colour_manual(values = c("#1f78b4","#e31a1c","#fdbf6f")) +
  scale_fill_manual(values = c("#1f78b4","#e31a1c","#fdbf6f")) +
  scale_y_discrete(limits=rev) +
  geom_ribbon(aes(xmin=mean-sd,xmax=mean+sd, fill = factor(Reading)), alpha=0.3) +
  theme_classic()
par_line2
ggsave('par_line.pdf', plot = par_line2, height = 4, width = 6, units = "in")


```
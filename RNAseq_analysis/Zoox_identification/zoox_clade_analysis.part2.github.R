# directory with the results from the part1
indir='/home/gospozha/haifa/cayman/rna/zoox/R/'
# design file
design_file='/home/gospozha/haifa/cayman/rna/zoox/design.all.csv'
design_factor='group'
# the order of factor levels
design_factor_levels=c("S", "D", "SS", "SD", "DD", "DS")
design_sample='sample'
# output directory
outidx='/home/gospozha/haifa/cayman/rna/zoox/R/'

# import necessary libraries
library(data.table)
library(plyr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(vegan)
library(pairwiseAdonis)

# read the design file 
t_design=read.csv(design_file,sep=";",stringsAsFactors=F)
t_design1 = t_design[,c(design_sample,design_factor)]
names(t_design1) = c('sample','group')
t_design1$group = factor(t_design1$group,levels=design_factor_levels)

# read the list of frequencies from the part1
speciesFreq <- readRDS(paste0(indir,"speciesFreq.RData"))

# transform the list to a dataframe
for(i in 1:length(speciesFreq)){
  if(i==1){
    speciesFreq1=speciesFreq[[i]]
  }else{
    speciesFreq1=merge(speciesFreq1,speciesFreq[[i]],by='species',all.x=TRUE,all.y=TRUE)
  }
}
# sort frequency table by the total counts in descending order
speciesFreq2=speciesFreq1[order(-apply(speciesFreq1[,2:ncol(speciesFreq1)],1,sum)),]

# save the frequency table to a csv
write.csv(speciesFreq2,file=paste0(outidx,"species-freq.csv"), row.names=FALSE)

# prepare the proportion table from the frequency table
prop_table <- speciesFreq2 %>%
  mutate(across(-species, ~ . / sum(.))) %>% 
  as.data.frame() 

# save the proportion table to a csv
write.csv(prop_table,file=paste0(outidx,"species-prop.csv"), row.names=FALSE)

# melt frequency table 
mlt1b=melt(speciesFreq2,id.vars=c('species'))
mlt2b=merge(mlt1b,t_design1,by.x='variable',by.y='sample',all.x=T,all.y=T)
# create a pdf canvas
pdf(paste0(outidx,"species-count-all.pdf"),width=10)
# ggplot2 heatmap
hm2 = ggplot(data = mlt2b, mapping = aes(x = variable, y = species, fill = value)) +
  geom_raster() + ggtitle('Count of mapped reads') +
  xlab("samples") +
  ylab("symbiont species") +  
  guides(fill=guide_legend(title="N of reads")) +
  theme(axis.text.x=element_text(angle = 45,size = 8, vjust=.8, hjust=0.8)) + 
  facet_grid(cols = vars(group),scale='free') +
  scale_fill_gradientn(colours = c('white','darkorange1','darkred'),oob = scales::squish)+
  scale_y_discrete(labels = c("Dtrenchii" = "Durusdinium trenchii", "symb" = "Breviolum minutum",
                              "syma" = "Symbiodinium clade A3", "Stri" = "Symbiodinium tridacnidorum", 
                              "SPilosum" = "Symbiodinium pilosum", "Snec" = "Symbiodinium necroappetens",
                              "Snat" = "Symbiodinium natans", "SMicUQSCI" = "Symbiodinium microadriaticum", 
                              "Slin" = "Symbiodinium linuacheae", "SKawaF" = "Fugacium kawagutii", 
                              "SGoreauC" = "Cladocopium goreau", "CladocPluteaC" = "Cladcopium C15"))
print(hm2)
dev.off()

# melt proportion table and draw a ggplot2 heatmap
mlt1=melt(prop_table,id.vars='species')
mlt2=merge(mlt1,t_design1,by.x='variable',by.y='sample',all.x=T,all.y=T)
pdf(paste0(outidx,"species-freq-all.pdf"),width=10)
hm1 = ggplot(data = mlt2, mapping = aes(x = variable, y = species, fill = value)) +
  geom_raster() + ggtitle('Proportion of mapped reads') +
  xlab("samples") +
  ylab("symbiont species") +  
  guides(fill=guide_legend(title="% of reads")) +
  theme(axis.text.x=element_text(angle = 45,size = 8, vjust=.8, hjust=0.8)) + 
  facet_grid(cols = vars(group),scale='free') +
  scale_fill_gradientn(colours = c('white','darkorange1','darkred'),limits=c(0,0.5), oob = scales::squish)+
  scale_y_discrete(labels = c("Dtrenchii" = "Durusdinium trenchii", "symb" = "Breviolum minutum",
                              "syma" = "Symbiodinium clade A3", "Stri" = "Symbiodinium tridacnidorum", 
                              "SPilosum" = "Symbiodinium pilosum", "Snec" = "Symbiodinium necroappetens",
                              "Snat" = "Symbiodinium natans", "SMicUQSCI" = "Symbiodinium microadriaticum", 
                              "Slin" = "Symbiodinium linuacheae", "SKawaF" = "Fugacium kawagutii", 
                              "SGoreauC" = "Cladocopium goreau", "CladocPluteaC" = "Cladcopium C15"))

print(hm1)
dev.off()


# the number of aligned reads
nofreads <- as.data.frame(colSums(speciesFreq2[,-1]))
write.csv(nofreads,file=paste0(outidx,"Nofalignedreads.csv"), row.names=TRUE)

# similarity indexes
design.ord <- t_design[order(t_design$sample),]
design.ord$depth = factor(design.ord$depth,levels=design_factor_levels)
design.ord$site = factor(design.ord$site,levels=c("CC", "MF"))
design.ord <- subset(design.ord, select = -c(sample,group) )

rownames(prop_table) <- prop_table$species
prop_table.bray <- t(prop_table[,-1])
bray_curtis <- vegdist(prop_table.bray, method = "bray")
pcoa_result <- cmdscale(bray_curtis, k = 2, eig = TRUE)

# Convert to a data frame for plotting
pcoa_df <- data.frame(Sample = rownames(prop_table.bray), 
                      Axis1 = pcoa_result$points[,1], 
                      Axis2 = pcoa_result$points[,2])

# Plot PCoA
ggplot(pcoa_df, aes(x = Axis1, y = Axis2, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -1) +
  theme_minimal() +
  labs(x = "PCoA Axis 1", y = "PCoA Axis 2", title = "PCoA - Bray-Curtis Dissimilarity")

rownames(t_design) <- t_design$sample
t_design$origin <- factor(t_design$origin)
t_design$depth <- factor(t_design$depth)

adonis2(bray_curtis ~ depth, data = t_design, permutations = 10000)
adonis2(bray_curtis ~ site * depth, data = t_design, permutations = 10000)
# both are not significant

# groups have equal dispersion, valid permanova results
dispersion <- betadisper(bray_curtis, t_design$depth)
permutest(dispersion)
dispersion <- betadisper(bray_curtis, t_design$site)
permutest(dispersion)

# pairwise adonis
pairwise.adonis(bray_curtis, t_design$depth, perm = 999)
# not significant

# compare group distances
# convert distance matrix to long dataframe
dist_df <- as.data.frame(as.table(as.matrix(bray_curtis)))
colnames(dist_df) <- c("Sample1", "Sample2", "Distance")

# add group info to each sample
dist_df <- dist_df %>%
  left_join(t_design %>% dplyr::select(Sample1 = sample, Depth1 = depth), by = "Sample1") %>%
  left_join(t_design %>% dplyr::select(Sample2 = sample, Depth2 = depth), by = "Sample2") %>%
  filter(Sample1 != Sample2) # remove self-comparisons

# tag comparisons as 'within' or 'between'
dist_df <- dist_df %>%
  mutate(Comparison = ifelse(Depth1 == Depth2, "Within", "Between"))

# compare the distances
wilcox.test(Distance ~ Comparison, data = dist_df)

simil <- ggplot(dist_df, aes(x = Comparison, y = Distance, fill = Comparison)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  annotate("text", x=2, y=0.8, label= "wilcox.test p.value=0.25", size=3) + 
  labs(title = "Within vs Between Group Dissimilarity", y = "Bray-Curtis Distance")
simil
ggsave("symb.bray.jpg", simil, width = 6, height = 6)  

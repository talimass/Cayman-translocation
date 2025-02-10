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
  scale_y_discrete(labels = c("symbd" = "Durusdinium sp.", "symb" = "Breviolum minutum",
                              "syma" = "Symbiodinium clade A3", "Stri" = "Symbiodinium tridacnidorum", 
                              "SPilosum" = "Symbiodinium pilosum", "Snec" = "Symbiodinium necroappetens",
                              "Snat" = "Symbiodinium natans", "SMicUQSCI" = "Symbiodinium microadriaticum UQ", 
                              "SMicUQCass" = "Symbiodinium microadriaticum UQ Cass", 
                              "SMicReefgen" = "Symbiodinium microadriaticum RG",
                              "Slin" = "Symbiodinium linuacheae", "SKawaF" = "Fugacium kawagutii", 
                              "SGoreauC" = "Cladocopium goreau", "CladocPluteaC" = "Cladcopium C15", "Cinf" = "Cladocopium infistulum"))
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
  scale_fill_gradientn(colours = c('white','darkorange1','darkred'),oob = scales::squish)+
  scale_y_discrete(labels = c("symbd" = "Durusdinium sp.", "symb" = "Breviolum minutum",
                              "syma" = "Symbiodinium clade A3", "Stri" = "Symbiodinium tridacnidorum", 
                              "SPilosum" = "Symbiodinium pilosum", "Snec" = "Symbiodinium necroappetens",
                              "Snat" = "Symbiodinium natans", "SMicUQSCI" = "Symbiodinium microadriaticum UQ", 
                              "SMicUQCass" = "Symbiodinium microadriaticum UQ Cass", 
                              "SMicReefgen" = "Symbiodinium microadriaticum RG",
                              "Slin" = "Symbiodinium linuacheae", "SKawaF" = "Fugacium kawagutii", 
                              "SGoreauC" = "Cladocopium goreau", "CladocPluteaC" = "Cladcopium C15", "Cinf" = "Cladocopium infistulum"))


print(hm1)
dev.off()


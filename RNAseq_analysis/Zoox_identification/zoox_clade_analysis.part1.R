# directory with DIAMOND results (*.m8 files)
indir='/lustre1/home/mass/eskalon/Porites/analysis/zoox/diamond.names.out'
files_pattern='(.*)(?=\\.m8$)'
# columns of *.m8 files
cols1=c('query' ,'subject',	'identity','alignment_length','mismatches','gap_openings','query_start','query_end','subject_start','subject_end','E_value','bit_score')
# output directory
outidx='/lustre1/home/mass/eskalon/Porites/analysis/zoox/R'

# import necessary libraries
library(data.table)
library(plyr)
library(ggplot2)
library(reshape2)

# make a list of *.m8 files
files0 = list.files(indir)
names0 = sapply(files0,function(x) regmatches(x,regexpr(files_pattern,x,perl=T))[1])
names1=names0[!is.na(names0)]
files1=files0[!is.na(names0)]

speciesFreq=list()

# for each file:
for(i in 1:length(files1)){
  # print the filename to stdout 
  cat(names1[i],"\n")
  f1=paste0(indir,"/",files1[i])
  # read the DIAMOND output file .m8
  blast1=fread(f1)
  # assign the column names
  names(blast1)=cols1
  # select only necessary columns
  blast1=blast1[,c('query' ,'subject','identity','alignment_length','E_value','bit_score')]
  # assign top priority to the top hit from DIAMOND
  blast1$i =1:nrow(blast1)
  top1=function(x){return(1:length(x))}
  blast1[, priority := top1(i),by=query]
  # keep only one top hit per query
  blast1 = blast1[priority ==1,]
  # match hit read to a certain species
  blast1[,species := regmatches(subject,regexpr('^[A-Za-z]+',subject,perl=TRUE))]
  cat("freq:\n")
  # compute and fill the list of species frequencies
  speciesFreq[[i]]=as.data.frame(table(blast1$species))
  # order by frequency
  speciesFreq[[i]]=speciesFreq[[i]][order(-speciesFreq[[i]]$Freq),]
  names(speciesFreq[[i]])=c('species',paste0(names1[i]))
}

# save the list of frequencies as an R object to use it in the next script
saveRDS(speciesFreq, file="speciesFreq.RData")
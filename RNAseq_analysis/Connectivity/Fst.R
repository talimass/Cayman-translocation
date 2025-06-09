###### Estimation of Fst statistic from transcriptome-derived SNPs ######

## Load libraries

library(gdsfmt)
library(SNPRelate)
library(hierfstat)
library(tidyr)
library(dplyr)
library(adegenet)
library(ade4)
#setwd("haifa/cayman/rna/connect/")
setwd("/lustre1/home/mass/eskalon/Porites/analysis/snp/Fst")

## Load the pruned SNPs

bed.fn <- "P_astreoides_vcf_pruned.bed"
fam.fn <- "P_astreoides_vcf_pruned.fam"
bim.fn <- "P_astreoides_vcf_pruned.bim"

## Create gds file
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "plink_pruned.gds")
snpgdsSummary("plink_pruned.gds")

# The file name: /home/gospozha/haifa/cayman/rna/connect/plink_pruned.gds 
# The total number of samples: 31 
# The total number of SNPs: 48473 
# SNP genotypes are stored in SNP-major mode (Sample X SNP).
# The number of valid samples: 31 

# The number of biallelic unique SNPs: 62 
## Open the gds file
genofile <- snpgdsOpen("plink_pruned.gds")

## Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

## Get population information
pop_code <- scan("pop.txt", what=character())

head(cbind(sample.id, pop_code)) #assumes the order of sample IDs is as the same as population codes

## Make a PCA to look at the genetic clusters
pca <- snpgdsPCA(genofile, autosome.only=FALSE)
#[1] Working space: 31 samples, 48,462 SNPs

# Make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

# Draw
pdf(paste0('PCA_prunedSNPs','.pdf'))
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))
dev.off()

# Look at the % variance explained by the components
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
#[1] 18.14  9.60  8.33  7.74  5.95  5.52

snpset <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
table(snpset)


x2<-2-snpgdsGetGeno(genofile) #snpgdsVCF2GDS stores the number of reference alleles, we want the the number of alternate alleles
#[1] Genotype matrix: 31 samples X 48473 SNPs

#add column with samples names 
vec <- sample.id

locus_transpose1 <- cbind(x2, sample = vec)  
locus_transpose2=as.data.frame(locus_transpose1)

#move last column with samples names to first
locus_transpose3 <- locus_transpose2 %>%
  select(sample, everything())

#add column with with groups as numbers, adult_meso is 1, planu_meso is 2 and so on
vec2 <- as.numeric(factor(pop_code))
# 5 2 2 5 5 2 4 1 1 4 1 4 1 1 4 4 4 1 6 3 6 2 3 3 3 2 3 6 6 6 5
# SD DD DD SD SD DD S D D S D S D D S S S D SS DS SS DD DS DS DS DD DS SS SS SS SD

# 5 SD
# 2 DD
# 4 S
# 1 D
# 6 SS
# 3 DS
#add column with depths
#vec3 <- as.numeric(factor(pop_code))

#vec2 <- c(1,1,1,2,2,1,2,1,2,1,2,1,2,1,2,2,1,1,2,2,1,1,2,1,1,1,1,1,2,2,2)
vec3 <- c(1,1,1,2,2,1,2,1,2,1,2,1,2,1,2,2,1,1,2,2,1,1,2,1,1,1,1,1,2,2,2)


# Combine as character labels (e.g. "5_1", "2_1", etc.)
combo_labels <- paste(vec2, vec3, sep = "_")

# Convert to numeric groupings
combo_numeric <- as.numeric(factor(combo_labels))

# View the result
combo_numeric

vec2 <- combo_numeric
vec3 <- combo_numeric
# 1 CC 2 MF
locus_transpose4 <- cbind(locus_transpose3, group = vec2) 
locus_transpose5 <- cbind(locus_transpose4, depth = vec3) 

#move last columns with depth and colony to 2nd and 3rd places

moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}

#move group column after sample column using the moveme function above
locus_transpose6=locus_transpose5[moveme(names(locus_transpose5), "group after sample")]

#move depth column after group column using the moveme function above
locus_transpose7=locus_transpose6[moveme(names(locus_transpose6), "depth after group")]

locus_transpose7$sample <- NULL #don't need this column

locus_transpose8 <- lapply(locus_transpose7,as.numeric) #pairwise.WCfst only works with numerics
locus_transpose9=as.data.frame(locus_transpose8)

Fst_between_groups = pairwise.WCfst(locus_transpose9[,-2],diploid=TRUE) #Fst between groups of samples

print(Fst_between_groups)
#1             2            3             4            5            6
#1           NA -0.0283448949 -0.040300063 -0.0026983709  0.020080942 -0.005871113
#2 -0.028344895            NA -0.046773579  0.0009893629  0.022110140 -0.007932272
#3 -0.040300063 -0.0467735787           NA -0.0131475093  0.025231677 -0.003746557
#4 -0.002698371  0.0009893629 -0.013147509            NA -0.001159314 -0.026099977
#5  0.020080942  0.0221101403  0.025231677 -0.0011593136           NA  0.004139822
#6 -0.005871113 -0.0079322715 -0.003746557 -0.0260999772  0.004139822           NA


# 5 SD
# 2 DD
# 4 S
# 1 D
# 6 SS
# 3 DS
## Permute Fst to get p values
#see the original script here https://lists.r-forge.r-project.org/pipermail/adegenet-forum/2011-February/000214.html

#1          2
#1         NA 0.07643047
#2 0.07643047         NA
# 1 CC 2 MF
# 
# 1             2             3             4             5             6             7             8            9           10            11
# 1             NA -0.0122809622 -3.651557e-02 -3.746043e-02  2.725063e-04  0.0022255332  0.0101968172  2.709894e-04 0.0241716291  0.003356072 -4.356488e-03
# 2  -0.0122809622            NA -1.221561e-02 -1.112792e-02  1.670942e-04  0.0006067508  0.0219450066  2.334843e-04 0.0265518087  0.004013930  5.052963e-04
# 3  -0.0365155662 -0.0122156117            NA -4.718307e-02  1.558172e-05  0.0066591771  0.0156390139  2.192639e-04 0.0242563520  0.004244788 -6.425526e-03
# 4  -0.0374604295 -0.0111279226 -4.718307e-02            NA -2.463263e-06  0.0026315918  0.0080411274  7.932744e-05 0.0270925743  0.004447467 -8.435873e-05
# 5   0.0002725063  0.0001670942  1.558172e-05 -2.463263e-06            NA -0.0004668597  0.0051105470 -9.917028e-05 0.0076146066  0.005764280  2.644712e-03
# 6   0.0022255332  0.0006067508  6.659177e-03  2.631592e-03 -4.668597e-04            NA  0.0624710021 -1.654070e-03 0.0583457829 -0.005809991  2.498063e-02
# 7   0.0101968172  0.0219450066  1.563901e-02  8.041127e-03  5.110547e-03  0.0624710021            NA  3.828090e-03 0.0001502513  0.026196465 -2.025846e-02
# 8   0.0002709894  0.0002334843  2.192639e-04  7.932744e-05 -9.917028e-05 -0.0016540698  0.0038280903            NA 0.0581830426  0.001917808  1.550344e-03
# 9   0.0241716291  0.0265518087  2.425635e-02  2.709257e-02  7.614607e-03  0.0583457829  0.0001502513  5.818304e-02           NA  0.031702959  8.164559e-03
# 10  0.0033560723  0.0040139302  4.244788e-03  4.447467e-03  5.764280e-03 -0.0058099913  0.0261964647  1.917808e-03 0.0317029586           NA  1.819923e-02
# 11 -0.0043564882  0.0005052963 -6.425526e-03 -8.435873e-05  2.644712e-03  0.0249806348 -0.0202584610  1.550344e-03 0.0081645586  0.018199230            NA


NBPERM <- 99 # this is the number of permutations used for the p-values
mat.perm <- lapply(1:NBPERM, function(i) pairwise.WCfst(locus_transpose9[,-2],diploid=TRUE))

#Fst_between_groups contains original Fst values, mat.perm is a list with NPERM matrices of permuted Fst values. 
#To get e.g. right-tail p-values, you can just count the proportion of mat.obs >= mat.perm; e.g. for the first pair of populations:
mean(c(Fst_between_groups[1,2] < na.omit(sapply(1:NBPERM, function(i) mat.perm[[i]][1,2])), TRUE))

#In the above command, "na.omit(sapply(1:NBPERM, function(i) mat.perm[[i]][1,2]))" is a vector of permuted values for this pair of populations across 
#all replicates; c(..., TRUE) is added because the observed value is always added to the permuted values (it is one of the possible permutations of groups).
#In practice, it is easier to convert the results as objects of the class randtest (class for Monte Carlo test in ade4):

Randtest <- as.randtest(na.omit(sapply(1:NBPERM, function(i) mat.perm[[i]][1,2])), Fst_between_groups[1,2], alter="greater")

#To have it done for all pairs of populations:
allTests <- list()
for(i in 1:(nrow(Fst_between_groups)-1)){
  for(j in 2:nrow(Fst_between_groups)){
    allTests[[paste(rownames(Fst_between_groups)[i],rownames(Fst_between_groups)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), Fst_between_groups[i,j], alter="greater")
  }
}  

print(allTests)

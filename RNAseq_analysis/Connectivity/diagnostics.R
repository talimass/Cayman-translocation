
library('gridExtra')
library('ggplot2')

# setting working directory 
setwd("/home/gospozha/haifa/cayman/rna/connect/")

VCFsnps <- read.csv('GVCFall_SNPs.table', header = T, na.strings=c("","NA"), sep = "\t") 
VCFindel <- read.csv('GVCFall_INDELs.table', header = T, na.strings=c("","NA"), sep = "\t")
dim(VCFsnps)
# 4142116     10
dim(VCFindel)
# 838000     10
VCF <- rbind(VCFsnps, VCFindel)
VCF$Variant <- factor(c(rep("SNPs", dim(VCFsnps)[1]), rep("Indels", dim(VCFindel)[1])))

snps <- '#A9E2E4'
indels <- '#F4CCCA'

DP <- ggplot(VCF, aes(x=DP, fill=Variant)) + geom_density(alpha=0.3) + 
  geom_vline(xintercept=c(10,6200)) 

QD <- ggplot(VCF, aes(x=QD, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=2, size=0.7)

FS <- ggplot(VCF, aes(x=FS, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(60, 200), size=0.7) + ylim(0,0.1)

MQ <- ggplot(VCF, aes(x=MQ, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=30, size=0.7)

MQRankSum <- ggplot(VCF, aes(x=MQRankSum, fill=Variant)) + geom_density(alpha=.3)

SOR <- ggplot(VCF, aes(x=SOR, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(4, 10), size=1, colour = c(snps,indels))

ReadPosRankSum <- ggplot(VCF, aes(x=ReadPosRankSum, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(-10,10,-20,20), size=1, colour = c(snps,snps,indels,indels)) + xlim(-30, 30)

pdf("Diag_plots.pdf", height=20, width=15)
theme_set(theme_gray(base_size = 18))
grid.arrange(QD, DP, FS, MQ, MQRankSum, SOR, ReadPosRankSum, nrow=4)
dev.off()

# Plots:
#	QUAL - variant confidence QUAL score (based on all reads)
#	QD - variant confidence standardized by depth.
#		- https://gatk.broadinstitute.org/hc/en-us/articles/360056968272-QualByDepth

# 	DP - *combined* depth per SNP across samples
#		- https://gatk.broadinstitute.org/hc/en-us/articles/360057441391-DepthPerSampleHC
#		- Depth of informative reads supporting variance (i.e. reads that show strong and clear evidence of the variant)

#	MQ - mapping quality of a SNP.
#		- https://gatk.broadinstitute.org/hc/en-us/articles/360057438331-RMSMappingQuality
#		- Overall mapping quality of reads supporting a variant call, averaged over all samples in a cohort.

#	FS - strand bias in support for REF vs ALT allele calls
#		- https://gatk.broadinstitute.org/hc/en-us/articles/360056968392-FisherStrand
#		- The output is a Phred-scaled p-value. The higher the output value, the more likely there is to be bias. More bias is indicative of false positive calls.

#	SOR - sequencing bias in which one DNA strand is favored over the other 
#		- https://gatk.broadinstitute.org/hc/en-us/articles/360056967852-AS-StrandOddsRatio
#		- Allele-specific strand bias estimated by the Symmetric Odds Ratio test.
 
#	MQRankSum - Rank sum test for mapping qualities of REF vs. ALT reads.
#		- https://gatk.broadinstitute.org/hc/en-us/articles/360057439151-MappingQualityRankSumTest
#		- Compares the mapping qualities of the reads supporting the reference allele with those supporting the alternate allele

#	ReadPosRankSum - do all the reads support a SNP call tend to be near the end of a read.
#		- https://gatk.broadinstitute.org/hc/en-us/articles/360057439971-AS-ReadPosRankSumTest
#		- Allele-specific Rank Sum Test for relative positioning of REF versus ALT allele within reads.
#		- Tests whether there is evidence of bias in the position of alleles within the reads that support them, between the reference and each alternate allele. 

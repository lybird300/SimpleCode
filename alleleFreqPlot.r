tagFreq <- read.table("tagSNPMAF.txt", header=TRUE, sep="\t")
hist(tagFreq$MinorAlleleFreq, breaks=10, main="Tag SNPs minor allele frequency distribution")
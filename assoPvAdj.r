args=(commandArgs(TRUE))
if (length (args) == 0) {
   print("No arguments supplied.")
   print("Please provide the following information (in order):")
   print("input=<absolute path to the input file>")
   print("output=<absolute path to the output file>")
   print("assoMethod=<the specific association testing method, either mach2dat or gcta>")
 } else {
   for (i in 1:length(args)) {
      eval (parse (text = args[[i]] ))
   }
 }

library(data.table)
if(assoMethod=="mach2dat"){
  assoResult<-fread(paste("awk '$1~/AFFECT/ && $2~/SNP/ {print $2,$12}' ", input))
  #assoResult<-fread("awk '$1~/AFFECT/ && $2~/SNP/ {print $2,$12}' asso_8308741_00908_Rep0.log")
}
if(assoMethod=="gcta"){
  assoResult<-fread(paste("awk 'NR>1{print $2,$9}' ", input))
}
assoResult<-as.data.frame(assoResult)
colnames(assoResult)=c("SNP","P")

## Bonferroni adjustment
bonferroni_p <- p.adjust(assoResult$P, method = "bonferroni")

## Holm adjustment
holm_p <- p.adjust(assoResult$P, method = "holm")

## Hochberg adjustment
hochberg_p <- p.adjust(assoResult$P, method = "hochberg")

## False discovery rate method q-value
fdr_q <- p.adjust(assoResult$P, method = "fdr")

adj.P <- data.frame(snp = assoResult$SNP,
                      nominal_p      = assoResult$P,
                      bonferroni_p = bonferroni_p,
                      holm_p       = holm_p,
                      hochberg_p   = hochberg_p,
                      fdr_q        = fdr_q)

#output snps at least one of whose adjusted p values is < 0.05
filtered.adj.P <- adj.P[which((adj.P$bonferroni_p<0.05) || (adj.P$holm_p<0.05) || (adj.P$hochberg_p<0.05) || (adj.P$fdr_q<0.05)), ]
if(nrow(filtered.adj.P)>0){
  write.table(filtered.adj.P, file=output)
}
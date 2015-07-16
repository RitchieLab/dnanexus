# usage: Rscript <this.R> <out_type>
args <- commandArgs(trailingOnly=TRUE)

# install.packages("/usr/share/SKAT_unlimited.tgz", repos=NULL, type="source")
library(SKAT)

bed.fn <- "input.bed"
bim.fn <- "input.bim"
fam.fn <- "input.fam"
cov.fn <- "input.covariate"
SKAT.setid <- "input-SetID.txt"
SSD.fn <- "temp.SSD"
INFO.fn <- "temp.INFO"
Generate_SSD_SetID(bed.fn, bim.fn, fam.fn, SKAT.setid, SSD.fn, INFO.fn)
FAM.data <- Read_Plink_FAM(fam.fn, Is.binary=FALSE)
FAM.cov <- Read_Plink_FAM_Cov(fam.fn, cov.fn, Is.binary=FALSE, cov_header=TRUE)

X = FAM.cov[,7:ncol(FAM.cov)]
obj.fn <- SKAT_Null_Model(FAM.data$Phenotype ~ as.matrix(X), out_type=args[1])
SSD.INFO <- Open_SSD(SSD.fn, INFO.fn)
output <- SKAT.SSD.All(SSD.INFO, obj.fn, kernel="linear")
output.mb <- SKAT.SSD.All(SSD.INFO, obj.fn, weights.beta=c(0.5,0.5))
output.default <- SKAT.SSD.All(SSD.INFO, obj.fn)

write.table(output$results, file="output-unweighted-uncorrected.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(output.mb$results, file="output-weighted-uncorrected.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(output.default$results, file="output-defaultweighted-uncorrected.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)

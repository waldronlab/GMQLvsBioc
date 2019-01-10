# "In TCGA data of BRCA patients, find the DNA somatic mutations 
# within the first 2000 bp outside of the genes that are both 
# expressed with FPKM > 3 and have at least a methylation in the same patient
# biospecimen, and extract these mutations of the top 5% patients
# with the highest number of such mutations."

library(curatedTCGAData)
system.time(mae <- curatedTCGAData("ACC", c("Mutation", "RNASeq2GeneNorm", "Methylation"), dry.run = FALSE)) 
## traditional methylation representation:
##    user  system elapsed 
## 171.516   3.109 181.086 
## HDF5 methylation representation (curatedTCGAData >= 1.5.6)
##  user  system elapsed 
## 6.223   0.223   9.607 

## create a data.frame for the methylation
meth <- assay(mae, "ACC_Methylation-20160128")
meth <- !is.na(meth) #NOT NA
meth <- meth * 1L #convert to integer mode
system.time(meth <- as.matrix(meth))
# number of non-missing per gene in each column
meth <- rowsum(meth, group=rowData(mae[["ACC_Methylation-20160128"]])$Gene_Symbol)
meth <- meth > 0 # logical "have at least a methylation"

## concatenate to mae. The only purpose of converting to 
## SummarizedExperiment here is to enable TCGAutils::symbolsToRanges(mae) below

mae <- c(mae, has.meth=SummarizedExperiment(meth))
rm(meth)

library(TCGAutils)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

mae <- TCGAutils::symbolsToRanges(mae)

## get ranges and add "within the first 2000 bp outside of the genes"
rr <- rowRanges(mae[["ACC_RNASeq2GeneNorm-20160128_ranged"]])
rr2000 <- flank(rr, width=2000, both=TRUE)

## summarize non-silent mutations per flanked range
nonsilent <- function(scores, ranges, qranges) any(scores != "Silent")
mutations <- RaggedExperiment::qreduceAssay(mae[["ACC_Mutation-20160128"]], rr2000, nonsilent, "Variant_Classification")
rownames(mutations) <- names(rr)
mutations[is.na(mutations)] <- 0
mae <- c(mae, mutations=SummarizedExperiment(mutations))
rm(mutations, rr, rr2000)

# select the elements for analysis
mae2 <- mae[, , c("ACC_RNASeq2GeneNorm-20160128_ranged", "has.meth_ranged", "mutations")]

saveRDS(mae2, file="ACC.rds")

mae3 <- intersectRows(mae2)
mae3 <- MatchedAssayExperiment(mae3)
keep <- assay(mae3, "ACC_RNASeq2GeneNorm-20160128_ranged") > 3 & 
  assay(mae3, "has.meth_ranged") &
  assay(mae3, "mutations")

summary(colSums(keep))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   5.000   8.000   8.961  11.000  69.000 

sort(colSums(keep), decreasing = TRUE)[1:round(ncol(keep)*0.05)]
# TCGA-PK-A5HB-01A-11R-A29S-07 TCGA-OR-A5LJ-01A-11R-A29S-07 
# 69                           29 
# TCGA-OR-A5J5-01A-11R-A29S-07 TCGA-OR-A5JA-01A-11R-A29S-07 
# 26                           22 

# > proc.time() with traditional methylation storage
#   almost all disk access (magnetic drive on my laptop), only ~4 seconds system time
#    user  system elapsed 
# 208.178   4.299 220.584 

# > proc.time() with HDF5 methylation storage
#   user  system elapsed 
# 35.139   2.641  48.287 

# "In TCGA data of BRCA patients, find the DNA somatic mutations 
# within the first 2000 bp outside of the genes that are both 
# expressed with FPKM > 3 and have at least a methylation in the same patient
# biospecimen, and extract these mutations of the top 5% patients
# with the highest number of such mutations."

library(curatedTCGAData)
system.time(mae <- curatedTCGAData("ACC", c("Mutation", "RNASeq2GeneNorm", "Methylation"), dry.run = FALSE)) #~3 minutes

## create a data.frame for the methylation
meth <- assay(mae, 1)
meth <- apply(meth, 2, as.numeric)  # soon unnecessary work-around
meth <- data.frame(meth, check.names = FALSE)
meth$symbols <- rowData(mae[["ACC_Methylation-20160128"]])$Gene_Symbol

## create a tibble that is TRUE if there is at least one 
## non-NA methylation value for the specimen
library(dplyr)
anyna <- function(x) any(is.na(x))
system.time(hasmeth <- dplyr::filter(meth, !is.na(symbols)) %>%
              dplyr::group_by(symbols) %>% 
              dplyr::summarize_all(anyna)) #129s user, 12.8s system, 142s elapsed

## convert to a matrix with symbols for rownames
hasmeth.matrix <- as.matrix(dplyr::select(hasmeth, -symbols))
rownames(hasmeth.matrix) <- hasmeth$symbols

## concatenate to mae
mae <- c(mae, has.meth=SummarizedExperiment(hasmeth.matrix))

library(TCGAutils)
mae <- TCGAutils::symbolsToRanges(mae)

## get ranges
rr <- rowRanges(mae[["ACC_RNASeq2GeneNorm-20160128_ranged"]])
rr2000 <- flank(rr, width=2000, both=TRUE)

## summarize non-silent mutations per flanked range
nonsilent <- function(scores, ranges, qranges) any(scores != "Silent")
mutations <- RaggedExperiment::qreduceAssay(mae[["ACC_Mutation-20160128"]], rr2000, nonsilent, "Variant_Classification")
rownames(mutations) <- names(rr)
mutations[is.na(mutations)] <- 0
mae <- c(mae, mutations=SummarizedExperiment(mutations))

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
# 0.000   3.000   5.000   6.727   8.000  55.000 

sort(colSums(keep), decreasing = TRUE)[1:5]
# TCGA-PK-A5HB-01A-11R-A29S-07 TCGA-OR-A5LJ-01A-11R-A29S-07 TCGA-OR-A5J5-01A-11R-A29S-07 
# 55                           26                           18 
# TCGA-OR-A5JA-01A-11R-A29S-07 TCGA-OR-A5JG-01A-11R-A29S-07 
# 18                           14 

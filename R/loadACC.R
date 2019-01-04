library(curatedTCGAData)
system.time(mae <- curatedTCGAData("ACC", c("Mutation", "RNASeq2GeneNorm", "Methylation"), dry.run = FALSE))
save(mae, file="~/ACC.rda")

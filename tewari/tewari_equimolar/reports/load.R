library(bcbioSmallRna)
meta = read.csv("metadata.csv", row.names = 1)
meta[["sample"]] = rownames(meta)
bcb = loadSmallRnaRun("final/2018-08-17_bcbio_eq-merged", 
                      interestingGroups = "lab",
                      colData = meta)
save(bcb, file = "data/bcb.rda")

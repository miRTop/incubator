library(tidyverse)

full_join(
    read_csv("gsm.txt"),
    read_tsv("SraRunTable.txt"), 
    by = c("gsm" = "Sample_Name")) %>%
    select(samplename = "Run", description = "name") %>% 
    write_csv("samples_bcbio_prepare.csv")

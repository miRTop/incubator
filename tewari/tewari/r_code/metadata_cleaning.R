setwd(here::here())

library(tidyverse)
library(rio)
library(janitor)

metadata = import("CrossU01_Metadata_Cleaned (1).xlsx") %>% clean_names() %>% 
    filter(grepl("QCPath.genome", file_type))

bcbio = list.files("tewari/bcbio/mirtop", pattern = "ready.gtf")

bcbio_meta = metadata %>% mutate(fixed_name = gsub("sample_", "", file_base_name),
                    fixed_name = gsub("_fastq", "", fixed_name)) %>% 
    filter(fixed_name  %in% gsub("-mirbase-ready.gtf", "", bcbio)) %>% 
    select(fixed_name, lab, lib_method_simple, lane, replicate) %>% 
    filter(lab != "Lab3",
           !grepl("1152v1", fixed_name),
           !grepl("plasma-[1234]-4N", fixed_name))

bcbio_meta %>% count(lab,lib_method_simple)

pilot_lab = paste0("Lab", c(1, 2, 4, 5))
pilot_lib = c("TruSeq", "NEBNext")

meta_pilot = bcbio_meta %>% 
    filter(lab  %in% pilot_lab,
           lib_method_simple  %in% pilot_lib) %>% 
    write_csv("meta_pilot.csv")



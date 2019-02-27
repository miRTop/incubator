setwd(here::here())
library(tidyverse)
library(ggplot2)
library(cowplot)
source("r_code/functions.R")
theme_set(theme_bw(base_size = 14))

meta_pilot = read_csv("meta_pilot.csv")

complete = read_tsv("tools/bcbio/mirtop/expression_counts.tsv.gz")

pilot = cbind(complete) %>%
    select(-Read, -iso_5p:-iso_snp_nt)  %>% 
    gather(sample, value, -UID, -miRNA, -Variant) %>% 
    mutate(sample = gsub("-mirbase-ready", "", sample)) %>%
    filter(sample  %in%  meta_pilot[["fixed_name"]]) %>%
    filter(value >= 1) %>%
    mutate(Variant = ifelse(is.na(Variant), "Reference", Variant)) %>%
    left_join(meta_pilot, by = c("sample" = "fixed_name"))


top_reference=function(r,v){
    sum(which(r==1) == which(v=="Reference"))
}

group_by(pilot, miRNA, sample) %>% 
    arrange(sample, miRNA, desc(value)) %>% 
    mutate(rank = 1:n()) %>% 
    summarise(total=n(), reference=top_reference(rank,Variant)) %>% 
    group_by(sample) %>% 
    summarise(mirnas=n(), isomirs=sum(total), reference=sum(reference)) %>% 
    left_join(meta_pilot, by = c("sample" = "fixed_name")) %>% 
    ggplot(aes(replicate, mirnas, fill="non-reference")) +
    geom_bar(stat="identity") +
    geom_bar(aes(replicate, reference, fill="reference"), stat="identity") +
    geom_text(aes(replicate, mirnas+50, label=isomirs), size=3) +
    facet_grid(lib_method_simple~lab) +
    ggsave("figures/mirna_stats/bcbio.pdf",
           width = 9, height = 7)
    

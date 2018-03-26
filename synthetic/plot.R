library(tidyverse)
library(ggplot2)
library(readr)

fn = c(synthetic = "synthetic_100/mirtop_stats.txt",
       bcbio = "bcbio/mirtop_stats.txt",
       mirge = "mirge/mirtop_stats.txt",
       srnabench = "sRNAbench/mirtop_stats.txt",
       isomirsea = "isomirsea/mirtop_stats.txt",
       prost = "prost/mirtop_stats.txt")

df = lapply(names(fn), function(x){
    read_csv(fn[x]) %>% 
        mutate(tool = x)
} ) %>% bind_rows()

ggplot(df %>%  filter(grepl("_sum", category)), aes(tool, counts)) +
    geom_bar(stat = "identity") + 
    facet_wrap(~category) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    ggsave("benchmark_sum.png")

ggplot(df %>%  filter(grepl("_count", category)), aes(tool, counts)) +
    geom_bar(stat = "identity") + 
    geom_text(aes(label = counts, y = counts + 5 )) +
    ylim(c(0, 50)) +
    facet_wrap(~category) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    ggsave("benchmark_counts.png")

read_tsv("all/summary.txt") %>% select(sample, idu, tag, same_mirna, iso_3p:iso_snp_central) %>% 
    gather("iso", "value", -sample, -idu, -tag, -same_mirna) %>%
    filter((tag == "D") | (tag == "E" & value == "FP") | (tag == "M" & value == "FN")) %>% 
    mutate(accuracy = value,
           accuracy = ifelse(tag == "E", "Extra", accuracy),
           accuracy = ifelse(tag == "M", "Miss", accuracy),
           accuracy = ifelse(tag == "D" & same_mirna != "Y", "CrossMap", accuracy)) %>% 
    distinct() %>% 
    filter(accuracy != "TN") %>% 
    ggplot(aes(x = sample, fill = accuracy)) + 
    geom_bar(position = "dodge")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    facet_wrap(~iso, scales = "free_y")  + 
    ggsave("benchmark_reference.png")

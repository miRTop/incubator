library(tidyverse)
library(ggplot2)
library(readr)

fn = c(synthetic = "synthetic/mirtop_stats.txt",
       bcbio = "tools/bcbio/mirtop_stats.txt",
       mirge = "tools/mirge/mirtop_stats.txt",
       srnabench = "tools/sRNAbench/mirtop_stats.txt",
       isomirsea = "tools/isomirsea/mirtop_stats.txt",
       prost = "tools/prost/mirtop_stats.txt")

df = lapply(names(fn), function(x){
    read_csv(fn[x]) %>% 
        mutate(tool = x)
} ) %>% bind_rows()

ggplot(df %>%  filter(grepl("_sum", category)), aes(tool, counts)) +
    geom_bar(stat = "identity") + 
    facet_wrap(~category) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    ggsave("plots/benchmark_sum.png")

ggplot(df %>%  filter(grepl("_count", category)), aes(tool, counts)) +
    geom_bar(stat = "identity") + 
    geom_text(aes(label = counts, y = counts + 5 )) +
    ylim(c(0, 50)) +
    facet_wrap(~category) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    ggsave("plots/benchmark_counts.png")

read_tsv("all/summary.txt") %>% select(sample, idu, tag, same_mirna, iso_3p:iso_snp_central) %>% 
    gather("iso", "value", -sample, -idu, -tag, -same_mirna) %>%
    filter((tag == "D") | (tag == "E" & value == "FP") | (tag == "M" & value == "FN")) %>% 
    mutate(tools = ifelse(grepl("sRNAbench", sample), "sRNAbench", "NA"),
           tools = ifelse(grepl("ready", sample), "bcbio", tools),
           tools = ifelse(grepl("0327", sample), "miRge", tools),
           tools = ifelse(grepl("tag", sample), "isomiRSEA", tools)) %>% 
    mutate(accuracy = value,
           accuracy = ifelse(tag == "E", "Extra", accuracy),
           accuracy = ifelse(tag == "M", "Miss", accuracy),
           accuracy = ifelse(tag == "D" & same_mirna != "Y", "CrossMap", accuracy)) %>% 
    distinct() %>% 
    filter(accuracy != "TN") %>% 
    write_csv("all/summary_parsed.txt")

read_csv("all/summary_parsed.txt") %>% 
    ggplot(aes(x = tools, fill = accuracy)) + 
    geom_bar(position = "dodge")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    facet_wrap(~iso, scales = "free_y")  + 
    scale_fill_brewer(palette = "Set1") +
    ggsave("plots/benchmark_reference.png")

read_csv("all/summary_parsed.txt") %>%
    mutate(iso = ifelse(iso == "iso_add", "iso_3p", iso),
           iso = ifelse(grepl("snp", iso), "iso_snp_all", iso)) %>% 
    ggplot(aes(x = tools, fill = accuracy)) + 
    geom_bar(position = "dodge")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    facet_wrap(~iso, scales = "free_y")  + 
    scale_fill_brewer(palette = "Set1") +
    ggsave("plots/benchmark_reference_simpler.png")

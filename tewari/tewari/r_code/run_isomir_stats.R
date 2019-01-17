setwd(here::here())
library(tidyverse)
library(ggplot2)
library(cowplot)
library(edgeR)
source("r_code/functions.R")
theme_set(theme_bw(base_size = 14))

tables = "tables/stats"

source("r_code/stats/bcbio_stats.R")
source("r_code/stats/mirge_stats.R")
source("r_code/stats/isomirsea_stats.R")
source("r_code/stats/srnabench_stats.R")

theme_set(theme_bw(base_size = 14))
stats = lapply(list.files(tables, full.names = TRUE), function(fn){
    tool = gsub("_stats_fixed_names.csv", "", basename(fn))
    read_csv(fn) %>%
        mutate(tool = tool)
}) %>% bind_rows()

stats %>%
    filter(grepl("_sum", category)) %>%
    mutate(category=ifelse(grepl("snp", category), "iso_snp", category),
           category = gsub("_sum", "", category)) %>%
    filter(category != "isomiR") %>%
    ggplot(aes(category, counts, color = tool, shape = as.factor(replicate))) +
    geom_point() +
    facet_grid(lab~lib_method_simple) +
    scale_color_brewer(palette = "Set2") +
    ggsave("figures/stats/summary_sum.png", width = 12, height = 6)


stats %>%
    filter(grepl("_count", category)) %>%
    mutate(category=ifelse(grepl("snp", category), "iso_snp", category),
           category = gsub("_count", "", category)) %>%
    ggplot(aes(category, counts, color = tool, shape = as.factor(replicate))) +
    geom_point() +
    facet_grid(lab~lib_method_simple) +
    scale_color_brewer(palette = "Set2")+
    ggsave("figures/stats/summary_counts.png", width = 12, height = 6)
